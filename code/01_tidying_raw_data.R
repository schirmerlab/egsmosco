suppressWarnings({
  library("reticulate")
  library("R.utils")
  library("ggplot2")
  library("SingleCellExperiment")
  library("scater")
  library("Seurat")
  library("Matrix")
  library("data.table")
  library("stringr")
})

mkdirs("output/01_tidy_data")

# Marques ####

# load metadata
marques_metadata <- read.csv(file = "data/marques_metada.csv")

# load matrix
marques_exp_mat <- readMM(file = "data/marques_exp_mat.mtx")

# load col and rownames
marques_exp_mat_colnames <- read.csv(file = "data/marques_exp_mat_colnames.csv")

marques_exp_mat_colnames <- marques_exp_mat_colnames$x
marques_exp_mat_rownames <- read.csv(file = "data/marques_exp_mat_rownames.csv")
marques_exp_mat_rownames <- marques_exp_mat_rownames$x

# asign col and rownames ( it is immportant because when saving sparse matrix row and colnames are delted)
colnames(marques_exp_mat) <- marques_exp_mat_colnames
rownames(marques_exp_mat) <- marques_exp_mat_rownames

# check whether ordering of cells in metada is the same as in matrix
all(marques_metadata$cell_barcode == colnames(marques_exp_mat))

# it is not, order them
indx <- match(colnames(marques_exp_mat), marques_metadata$cell_barcode)
marques_metadata <- marques_metadata[indx,]

# check now
all(marques_metadata$cell_barcode == colnames(marques_exp_mat))

#add rownames in metadata (Seurat requires that)
rownames(marques_metadata) <- marques_metadata$cell_barcode

marques <- CreateSeuratObject(counts = marques_exp_mat, meta.data = marques_metadata)
marques <- SCTransform(marques)
marques <- RunPCA(marques)
ElbowPlot(marques)
marques <- RunUMAP(marques, dims = 1:18) # chosen based on get_sig_pc (not shown)
marques <- FindNeighbors(marques, dims = 1:18)
marques <- FindClusters(marques)

Idents(marques) <- marques$cell_subtypes
UMAPPlot(marques)
marques <- RenameIdents(marques, "PPR" = "VLMC")
UMAPPlot(marques)
marques$cell_types_short <- Idents(marques)
marques$cell_types_short <- factor(marques$cell_types_short,
                                   levels = sort(levels(marques$cell_types_short)))
Idents(marques) <- marques$cell_types_short
DimPlot(marques)
saveRDS(marques, file = "output/01_tidy_data/marques.rds")


# Floridia ####
load("data/floriddia.Rdata")
Idents(oligos.integrated) <- oligos.integrated$confidentclusters
UMAPPlot(oligos.integrated)
floriddia <- oligos.integrated

floriddia[["cell_types_short"]] <- Idents(floriddia)
floriddia$cell_types_short <- factor(floriddia$cell_types_short,
                                     levels = sort(levels(floriddia$cell_types_short)))
UMAPPlot(floriddia)
saveRDS(floriddia, file = "output/01_tidy_data/floriddia.rds")


# Schirmer ####
# Load normalized expression matrix. It is in tsv file and it is huge read it with fread

schirmer_norm_mtx <-
  fread(file = "data/schirmer_norm_data/schirmer_norm_mat.tsv",
        sep = "\t",
        data.table = F)

head(rownames(schirmer_norm_mtx))
head(colnames(schirmer_norm_mtx))

# gene names which should be rownames of the matrix are stored in a separate column named "gene" within the matrix

head(schirmer_norm_mtx$gene)
# in addition there are both ensable and gene symbols in the same column separated by "|".
# Wrangle that and obtain gene symbols only
genes_combined <- schirmer_norm_mtx$gene
gene_symbols <- str_extract(genes_combined, "(?<=\\|).*")

# now that i have gene symbols, i will delete column "gene" in the data frame, and will set rownames with gene_symbols that i extracted

# delete the "gene" column
schirmer_norm_mtx$gene <- NULL

# check duplication
sum(duplicated(gene_symbols)) #it is duplicated...

# since genes are duplicated i will name each duplicated instance with separate index and asign that to matrix
rownames(schirmer_norm_mtx) <- make.names(gene_symbols, unique = T)

#convert it to matrix
schirmer_counts <- as.matrix(schirmer_norm_mtx)

#convert it to sparse matrix
schirmer_counts <- Matrix(schirmer_counts, sparse = T)

# load metadta now
schirmer_metadata <- read.delim("data/schirmer_raw_data/meta.txt")

#check whether cells in metadata matches cells in matrix
all(schirmer_metadata$cell == colnames(schirmer_counts))

# it matches (if it wouldnt you would need to wrangle them by indexing). Now set rownames of metadata with cell barcodes. Seting rownames in meatadat is required step by seurat
rownames(schirmer_metadata) <- schirmer_metadata$cell

# Create seurat object now
schirmer <- CreateSeuratObject(counts = schirmer_counts, meta.data = schirmer_metadata)

# Now I have a seurat object in which normalized counts are set in the count slot. I will put the same matrix in the data slot so that they are exactly at the place they should be :D
# this step is required because you cant create Seurat object directly with data slot (at least i dont know the way). So first you need to create the object with any matrix which will be stored in the count slot, and then once you have object you can store various matrix in various slots

schirmer <- SetAssayData(schirmer, "data", schirmer_counts)

# I will once again add metadata (not necesary since i created object using metadata. However, I will just using their AddMetadata function which is default way for stroing metadata)
schirmer <- AddMetaData(object = schirmer, metadata = schirmer_metadata)


# now I need to set embedings whic h are in the raw metadata
embedings <-
  schirmer_metadata %>%
  dplyr::select(tsne1, tsne2) %>%
  dplyr::rename(tSNE_1 = tsne1,
                tSNE_2 = tsne2) %>%
  as.matrix()

schirmer[["tsne"]] <- CreateDimReducObject(embedings, key = "tSNE_")

DimPlot(schirmer, reduction = "tsne", group.by = "cell_type")

TSNEPlot(schirmer)
Idents(schirmer) <- schirmer$cell_type
TSNEPlot(schirmer, cols = DiscretePalette(22), label =T)

new_cluster_ids <- c('EN', 'EN', 'IN', 'EN', 'EN', 'IN', 'IN', 'EC', 'OPC', 'EN', 'IN', 'EN', 'OL', 'AS', 'OL', 'SC', 'OL', 'MIX', 'MG', 'PHA', 'LYM', 'LYM')
names(new_cluster_ids) <- levels(schirmer)

schirmer <- RenameIdents(schirmer, new_cluster_ids)
TSNEPlot(schirmer)
schirmer[["cell_subtypes_short"]] <- Idents(schirmer)
saveRDS(object = schirmer, file = "output/01_tidy_data/schirmer.rds")
