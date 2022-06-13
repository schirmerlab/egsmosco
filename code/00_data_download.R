library(R.utils)
library(GEOquery) # for getting GEO data
library(tidyverse)
library(Matrix)
library(openxlsx)

# compiled 22.06.2020

# Blum et al. ####

# data from: Blum et al. Single-cell transcriptomic analysis of the adult mouse spinal cord reveals molecular diversity of autonomic and skeletal motor neurons. Nat Neurosci. 2021 Feb 15. doi: 10.1038/s41593-020-00795-0. Epub ahead of print. PMID: 33589834.

# downalod
download.file(url = "http://spinalcordatlas.org/assets/downloads/allexpsctl.h5ad",
              destfile = "data/blum_full.h5ad",
              method = "curl")

# Floriddia et al. ####
# data from: Floriddia et al. Distinct oligodendrocyte populations have spatial preference and different responses to spinal cord injury. Nat Commun. 2020 Nov 17;11(1):5860. doi: 10.1038/s41467-020-19453-x. PMID: 33203872; PMCID: PMC7673029.

#download
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE128525&format=file&file=GSE128525%5FSeuratObjSpCrdGEO%2ERdata%2Egz",
              destfile = "data/floriddia.Rdata.gz",
              method = "curl")

#unzip
gunzip("data/floriddia.Rdata.gz")


# Marques et al ####
# data from: Marques et al. Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system. Science. 2016 Jun 10;352(6291):1326-1329. doi: 10.1126/science.aaf6463. PMID: 27284195; PMCID: PMC5221728

# Get full Geo object
marques <- getGEO(GEO = "GSE75330")

# obatin metadata
marques_metada <- pData(phenoData(marques[[1]]))
head(marques_metada)
rownames(marques_metada) <- NULL

# subset metadata to include only rlevavnt info
marques_metada <-
  marques_metada %>%
  select(title,
         source_name_ch1,
         characteristics_ch1,
         characteristics_ch1.1,
         characteristics_ch1.2,
         characteristics_ch1.3,
         characteristics_ch1.4) %>%
  rename(cell_barcode = title,
         cns_region = source_name_ch1,
         strain = characteristics_ch1,
         age = characteristics_ch1.1,
         sex = characteristics_ch1.2,
         treatment = characteristics_ch1.3,
         cell_subtypes = characteristics_ch1.4)

head(marques_metada)

# tidy metadata
marques_metada$strain <-
  str_extract(string = marques_metada$strain,
              pattern = "(?<=strain: )[A-Za-z0-9]*") %>%
  na_if("") # to replace empty characters with NA

marques_metada$age <-
  str_extract(string = marques_metada$age,
              pattern = "(?<=age: )[A-Za-z0-9]*") %>%
  na_if("") # to replace empty characters with NA

marques_metada$sex <-
  str_extract(string = marques_metada$sex,
              pattern = "(?<=Sex: )[A-Za-z0-9]*") %>%
  na_if("") # to replace empty characters with NA

marques_metada$treatment <-
  str_extract(string = marques_metada$treatment,
              pattern = "(?<=treatment: )[A-Za-z0-9]*") %>%
  na_if("") # to replace empty characters with NA

marques_metada$cell_subtypes <-
  str_extract(string = marques_metada$cell_subtypes,
              pattern = "(?<=type: )[A-Za-z0-9]*") %>%
  na_if("") # to replace empty characters with NA

head(marques_metada)

# Obtain expression matrix
getGEOSuppFiles(GEO = "GSE75330", baseDir = "data/", makeDirectory = F)

# unzip it
gunzip(filename = "data/GSE75330_Marques_et_al_mol_counts2.tab.gz")

# read it
marques_exp_mat <- read.delim(file = "data/GSE75330_Marques_et_al_mol_counts2.tab",
                              sep = "\t",
                              header = T)

# set rownames
marques_exp_mat <- marques_exp_mat %>% column_to_rownames(var = "cellid")

# convert to matrix format and sparse format
marques_exp_mat <- as.matrix(marques_exp_mat)
marques_exp_mat <- Matrix(marques_exp_mat, sparse = T)

# check whether cells in matrix matches with cells in metadata
colnames(marques_exp_mat) %in%  marques_metada$cell_barcode  # Nothing matches ...

# It is not matching because in the matrix cells.ids are separated by "." while in the metadata they are separated with "-".

# replacing "-" with "." in metadata.
marques_metada$cell_barcode <- marques_metada$cell_barcode %>% str_replace_all(pattern = "-", replacement = "\\.")

# Check matching now
all(colnames(marques_exp_mat) %in%  marques_metada$cell_barcode) # not full match
sum(!colnames(marques_exp_mat) %in%  marques_metada$cell_barcode) # 16 from the matrix are not present in the metadata

# check which cells from matrix are not present in metadata
mtx_indx <- which(FALSE == colnames(marques_exp_mat) %in%  marques_metada$cell_barcode)

# cells present in mtx but not in metadata
mtx_colnmaes_not_in_metadata <- colnames(marques_exp_mat)[mtx_indx]

# check which cells from  metadata are not present in matrix
df_indx <- which(FALSE == marques_metada$cell_barcode %in% colnames(marques_exp_mat))

# cells present in metadata but not in matrix
metadata_cells_not_in_matrix <- marques_metada$cell_barcode[mtx_indx]

supplementary_xy <- tibble("Present in matrix but not metadata" = mtx_colnmaes_not_in_metadata,
                           "Present in metadata but not in matrix" = metadata_cells_not_in_matrix)

write.xlsx(x = supplementary_xy, file = "data/Marques_discrepancy.xlsx")

# in order to have matching matrix and metadata, intersection of cells will be taken and subset with it both matrix and metadata
shared_barcodes <- intersect(colnames(marques_exp_mat), marques_metada$cell_barcode)

marques_metada <- marques_metada %>% filter(cell_barcode %in% shared_barcodes)
marques_exp_mat <- marques_exp_mat[,shared_barcodes]

write_csv(x = marques_metada, file = "data/marques_metada.csv")
writeMM(obj = marques_exp_mat, file = "data/marques_exp_mat.mtx")

mtx_cols <- colnames(marques_exp_mat)
mtx_rows <- rownames(marques_exp_mat)

write.csv(mtx_cols, file = "data/marques_exp_mat_colnames.csv")
write.csv(mtx_rows, file = "data/marques_exp_mat_rownames.csv")

# Schirmer et al ####

# To reproduce Schirmer original object I need nomralized expression matrix and metadata which includes tSNE coordinates and cell_annotation
# tSNE coordinates are included in metadata file in rawMatrix.zip folder. I will donwalod rawMatrix folder just becuase I need that metadata file with tSNE coordinates. For everything else I will use normalized expression matrix.

# Download rawMatrix folder where is metadata with tSNE  coordinates
download.file(url = "https://cells.ucsc.edu/ms/rawMatrix.zip",
              destfile = "data/schirmer_raw_data.zip",
              method = "curl")

unzip(zipfile = "data/schirmer_raw_data.zip",
      exdir = "data/schirmer_raw_data")

# Download normalized expression matrix
mkdirs("data/schirmer_norm_data")
download.file(url = "https://cells.ucsc.edu/ms/exprMatrix.tsv.gz",
              destfile = "data/schirmer_norm_data/schirmer_norm_mat.tsv.gz",
              method = "curl")

gunzip("data/schirmer_norm_data/schirmer_norm_mat.tsv.gz")
