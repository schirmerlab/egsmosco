library(Seurat)
library(tidyverse)
library(patchwork)
library(kableExtra)

if (!file.exists("output/03_GLUL_dge")) {
  dir.create("output/03_GLUL_dge")
}


atxm <- theme(plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
              legend.title = element_text(size = 18, hjust = 0.5),
              legend.title.align = 0.5,
              legend.text = element_text(size = 21),
              plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
              legend.margin = margin(0,0,0,0),
              axis.text = element_text(size = 24),
              axis.title = element_text(size = 30),
              axis.ticks.length = unit(0.12, "cm"),
              axis.line = element_line(size = 0.7))


schirmer <- readRDS("output/01_tidy_data/schirmer.rds")
floriddia <- readRDS("output/01_tidy_data/floriddia.rds")

Idents(schirmer) <- schirmer$cell_subtypes_short
degs_schirmer <- VlnPlot(schirmer,
                         idents = c("OL", "AS"),
                         split.by = "diagnosis",
                         features = "GLUL",
                         cols = c("#00BFC4", "#F8766D")) +
  ylim(0,7.2) +
  theme(plot.title = element_text(face = c("italic")))
degs_schirmer

schirmer$cell_condition <- paste(schirmer$cell_subtypes_short, schirmer$diagnosis, sep = "_")
Idents(schirmer) <- schirmer$cell_condition

glul_ol <- FindMarkers(object = schirmer, ident.1 = "OL_Control", ident.2 = "OL_MS", features = "GLUL")
glul_ol <-
  glul_ol %>%
  add_column("Cell type" = "Oligodendrocytes",
             Comparison = "Control vs MS",
             Gene = "GLUL") %>%
  select(Gene, `Cell type`, Comparison, everything())
glul_as <- FindMarkers(object = schirmer, ident.1 = "AS_Control", ident.2 = "AS_MS", features = "GLUL")
glul_as <- glul_as %>%
  add_column("Cell type" = "Astrocytes",
             Comparison = "Control vs MS",
             Gene = "GLUL") %>%
  select(Gene, `Cell type`, Comparison, everything())

supplementary_table_xyz <- bind_rows(glul_ol, glul_as)
rownames(supplementary_table_xyz) <- NULL


supplementary_table_xyz %>%
  kbl() %>%
  kable_styling()
openxlsx::write.xlsx(x = supplementary_table_xyz,
                     file = "output/03_GLUL_dge/supplementary_table_xyz.xlsx",
                     overwrite = T)


p1 <- degs_schirmer +
  labs(x = "Cell types") +
  theme(plot.title = element_text(face = "italic")) +
  geom_line(data = tibble(x=c(0.8,1.2), y=c(6.5,6.5)), aes(x=x,y=y), inherit.aes = F) +
  geom_line(data = tibble(x=c(1.8,2.2), y=c(6.8,6.8)), aes(x=x,y=y), inherit.aes = F) +
  geom_text(data = tibble(x=c(1), y=c(6.7)), aes(x=x,y=y, label = "***"),size = 8, inherit.aes = F) +
  geom_text(data = tibble(x=c(2), y=c(7.2)), aes(x=x,y=y, label = "n.s."),size = 8, inherit.aes = F) +
  atxm +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = unit(c(0,0,0,2), c("cm","cm","cm","cm")))
p1

Idents(schirmer) <- schirmer$cell_subtypes_short
p2 <- DimPlot(schirmer,
              group.by = "diagnosis",
              pt.size = 0.5,
              ) +
  scale_shape_manual(values=seq(0,11)) +
  ggtitle("Schirmer et al, 2019") +
  atxm
p2

Idents(floriddia) <- floriddia$cell_types_short
degs_floriddia <- VlnPlot(floriddia,
                          idents = c("MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "MOL6"),
                          split.by = "orig.ident",
                          features = "Glul",
                          cols = c("#00BFC4", "#F8766D", "#00BA38")) +
  atxm +
  ylim(0,6)
degs_floriddia

floriddia$cell_condition <- paste(floriddia$cell_types_short, floriddia$orig.ident, sep = "_")
Idents(floriddia) <- floriddia$cell_condition

mol6_ctrl_vs_is <- FindMarkers(object = floriddia,
                                features = "Glul",
                                ident.1 = "MOL6_CTRL",
                                ident.2 = "MOL6_IS",
                                logfc.threshold = 0.1)

mol6_ctrl_vs_wd <- FindMarkers(object = floriddia,
                               features = "Glul",
                               ident.1 = "MOL6_CTRL",
                               ident.2 = "MOL6_WD",
                               logfc.threshold = 0.1)


p3 <- degs_floriddia +
  labs(x = "Cell types") +
  theme(plot.title = element_text(face = "italic")) +
  geom_line(data = tibble(x=c(5.7,6), y=c(5.3,5.3)), aes(x=x,y=y), inherit.aes = F) +
  geom_text(data = tibble(x=c(5.85), y=c(5.4)), aes(x=x,y=y, label = "***"), inherit.aes = F)
p3

p4 <- UMAPPlot(floriddia, group.by = "cell_types_short") + ggtitle("Floriddia et al, 2020")
p4

Idents(floriddia) <- floriddia$cell_types_short
floriddia <- RenameIdents(floriddia, "MOL1" = "MOL", "MOL2" = "MOL", "MOL3" = "MOL", "MOL4" = "MOL", "MOL5" = "MOL", "MOL6" = "MOL")
floriddia$pooled_oligos <- Idents(floriddia)

Idents(floriddia) <- floriddia$pooled_oligos
p5 <- DimPlot(floriddia) +
  ggtitle("Floriddia et al, 2020") +
  atxm
p5

degs_floriddia <- VlnPlot(floriddia,
                          idents = "MOL",
                          split.by = "orig.ident",
                          features = "Glul",
                          cols = c("#00BFC4", "#F8766D", "#00BA38")) +
  ylim(0,6) +
  atxm +
  theme(plot.title = element_text(face = "italic"))
degs_floriddia

floriddia$cell_condition <- paste(floriddia$pooled_oligos, floriddia$orig.ident, sep = "_")
Idents(floriddia) <- floriddia$cell_condition

mol_ctrl_vs_is <- FindMarkers(object = floriddia,
                               features = "Glul",
                               ident.1 = "MOL_CTRL",
                               ident.2 = "MOL_IS",
                               logfc.threshold = 0.1)

p6 <- degs_floriddia +
  labs(x = "Cell types") +
  theme(plot.title = element_text(face = "italic")) +
  geom_line(data = tibble(x=c(0.7,1), y=c(5.3,5.3)), aes(x=x,y=y), inherit.aes = F) +
  geom_text(data = tibble(x=c(0.85), y=c(5.4)), aes(x=x,y=y, label = "***"), size = 8, inherit.aes = F)
p6

design <- "
AB
CC
"

fig <-
  p2 +
  p1 +
  plot_spacer() +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = design,
              heights = unit(c(15,1), c("cm", "null"))) &
  theme(plot.tag = element_text(size = 56))


fig
ggsave(filename = "output/03_GLUL_dge/glul_dge_analysis.pdf",
       plot = fig,
       device = "pdf",
       width = 21*3,
       height = 29.7*3,
       units = "cm")


ggsave(filename = "output/03_GLUL_dge/glul_dge_analysis.png",
       plot = fig,
       width = 21*3,
       height = 6.8*3,
       units = "cm")
