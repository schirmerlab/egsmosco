suppressWarnings({
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(pals)
  library(R.utils)
  library(magick)
})

mkdirs("output/02_GLUL_expression")
mkdirs("images")

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



object_paths <- Sys.glob(filePath("output/01_tidy_data/", "*"))
object_names <- tools::file_path_sans_ext(basename(object_paths))

object_list <-
  object_paths %>%
  map(readRDS) %>%
  setNames(object_names)

list2env(object_list, envir = .GlobalEnv)
rm(object_list)

p1 <- DimPlot(marques, cols = glasbey(13), pt.size = 0.7) +
  ggtitle("Dataset from Marques et al., 2016") +
  atxm +
  theme(legend.position = c(0.87,0.7))
p1

p2 <- VlnPlot(marques, c("Glul", "Mobp", "Pdgfra"), cols = glasbey(13),ncol = 1, pt.size = 0.000001)
p2 <-
  p2 &
  atxm &
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = c("italic")))
p2[[2]] <- p2[[2]] + theme(axis.title.y = element_text(size = 30))

p2[[3]] <- p2[[3]] + theme(axis.text.x = element_text(angle = 45),
                           plot.margin = unit(c(0.3,0.3,0.3,2), "cm"))


design <-  "
AB
CC
"

fig1 <-
  p1 +
  p2 +
  plot_spacer() +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.5,1),
              heights = unit(c(16,1), c("cm", "null"))) &
  theme(plot.tag = element_text(size = 56))

fig1[[2]][[2]] <- fig1[[2]][[2]] + theme(plot.tag = element_blank())
fig1[[2]][[3]] <- fig1[[2]][[3]] + theme(plot.tag = element_blank())

fig1

ggsave(filename = "output/02_GLUL_expression/marques.pdf",
       plot = fig1,
       width = 210*3,
       height = 297*3,
       units = "mm")

ggsave(filename = "output/02_GLUL_expression/marques.png",
       plot = fig1,
       width = 210*3,
       height = 75*3,
       units = "mm")

img = image_read("output/02_GLUL_expression/marques.png")
img <- image_scale(image_scale(img,"33%"))
image_write(img, path = "images/marques.png", format = "png")

p3 <- DimPlot(blum, cols = glasbey(13), pt.size = 0.7) +
  ggtitle("Dataset from Blum et al., 2021") +
  atxm +
  theme(legend.position = c(0.8,0.2))
p3

p4 <-
  FeaturePlot(blum, c("Glul", "Mobp", "Aqp4", "Flt1","Cx3cr1", "Rbfox3")) &
  atxm &
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey20", size = 1),
        axis.line = element_blank(),
        plot.title = element_text(face = c("italic")))
p4[[1]] <- p4[[1]] + theme(plot.margin = unit(c(0.3,0.3,0.3,2), "cm"))

design <-  "
AB
CC
"

fig2 <-
  p3 +
  p4 +
  plot_spacer() +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1,1),
              heights = unit(c(22,1), c("cm", "null"))) &
  theme(plot.tag = element_text(size = 56))

for (i in 2:6) {
  fig2[[2]][[i]] <- fig2[[2]][[i]] + theme(plot.tag = element_blank()) # remove tags from plots 2:6
}


fig2
ggsave(filename = "output/02_GLUL_expression/blum.pdf",
       plot = fig2,
       width = 210*3,
       height = 297*3,
       units = "mm")

ggsave(filename = "output/02_GLUL_expression/blum.png",
       plot = fig2,
       width = 210*3,
       height = 90*3,
       units = "mm")

img = image_read("output/02_GLUL_expression/blum.png")
img <- image_scale(image_scale(img,"33%"))
image_write(img, path = "images/blum.png", format = "png")

p5 <- DimPlot(floriddia, cols = glasbey(14), pt.size = 0.68) +
  ggtitle("Dataset from Floriddia et al. 2020") +
  atxm +
  theme(legend.position = c(0.87,0.72))
p5
p6 <- VlnPlot(floriddia, c("Glul", "Mobp", "Pdgfra"), cols = glasbey(14),ncol = 1, pt.size = 0.000001)
p6
p6 <-
  p6 &
  atxm &
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = c("italic")))
p6[[2]] <- p6[[2]] + theme(axis.title.y = element_text(size = 30))

p6[[3]] <- p6[[3]] + theme(axis.text.x = element_text(angle = 45),
                           plot.margin = unit(c(0.3,0.3,0.3,2), "cm"))

design <-  "
AB
CC
"

fig3 <-
  p5 +
  p6 +
  plot_spacer() +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.5,1),
              heights = unit(c(16,1), c("cm", "null"))) &
  theme(plot.tag = element_text(size = 56))

fig3[[2]][[2]] <- fig3[[2]][[2]] + theme(plot.tag = element_blank())
fig3[[2]][[3]] <- fig3[[2]][[3]] + theme(plot.tag = element_blank())


fig3
ggsave(filename = "output/02_GLUL_expression/floriddia.pdf",
       plot = fig3,
       width = 210*3,
       height = 297*3,
       units = "mm")

ggsave(filename = "output/02_GLUL_expression/floriddia.png",
       plot = fig3,
       width = 210*3,
       height = 75*3,
       units = "mm")

img = image_read("output/02_GLUL_expression/floriddia.png")
img <- image_scale(image_scale(img,"33%"))
image_write(img, path = "images/floriddia.png", format = "png")

p7 <- DimPlot(schirmer, cols = glasbey(13), pt.size = 0.7) +
  ggtitle("Dataset from Schirmer et al., 2019") +
  atxm +
  theme(legend.position = c(0.91,0.19))
p7

p8 <-
  FeaturePlot(schirmer, c("GLUL", "MOBP", "AQP4", "DOCK8","GAD1", "SV2B")) &
  atxm &
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey20", size = 1),
        axis.line = element_blank(),
        plot.title = element_text(face = c("italic")))
p8[[1]] <- p8[[1]] + theme(plot.margin = unit(c(0.3,0.3,0.3,2), "cm"))

design <-  "
AB
CC
"

fig4 <-
  p7 +
  p8 +
  plot_spacer() +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1,1),
              heights = unit(c(22,1), c("cm", "null"))) &
  theme(plot.tag = element_text(size = 56))

for (i in 2:6) {
  fig4[[2]][[i]] <- fig4[[2]][[i]] + theme(plot.tag = element_blank()) # remove tags from plots 2:6
}


fig4
ggsave(filename = "output/02_GLUL_expression/schirmer.pdf",
       plot = fig4,
       width = 210*3,
       height = 297*3,
       units = "mm")

ggsave(filename = "output/02_GLUL_expression/schirmer.png",
       plot = fig4,
       width = 210*3,
       height = 92*3,
       units = "mm")

img = image_read("output/02_GLUL_expression/schirmer.png")
img <- image_scale(image_scale(img,"33%"))
image_write(img, path = "images/schirmer.png", format = "png")


