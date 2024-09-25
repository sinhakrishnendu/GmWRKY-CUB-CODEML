library(tidyverse)
library(ggplot2)
library(pheatmap)
library(scales)
library(viridis)
library(readxl)
###
rscu <- read_excel('RSCU_GmWRKY179.xlsx')
norm <- as.data.frame(lapply(rscu[-c(1)],rescale))
mat <- as.matrix(rscu[-c(1)])
pheatmap(as.matrix(norm), 
         cutree_cols = 5, 
         cluster_rows = FALSE, 
         color = colorRampPalette(rocket(100))(100),
         filename = 'RSCUHeatMap1.1.tiff',
         dpi = 300)
pheatmap(as.matrix(norm), 
         cutree_cols = 5, 
         cluster_rows = FALSE, 
         color = colorRampPalette(c("yellow", "darkblue"))(100),
         filename = 'RSCUHeatMap.tiff',
         dpi = 300)
grmat_r <- as.matrix(rscu[-c(1)])
pheatmap(mat, cutree_cols = 8, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
# Example with clustering
pheatmap(mat, scale = "row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")



rscut <- read_excel('RSCU179.xlsx')
normt <- as.data.frame(lapply(rscut[-c(1)],rescale))
matt <- as.matrix(rscut[-c(1)])
pheatmap(as.matrix(normt), 
         cutree_cols = 5, 
         cluster_rows = FALSE, 
         color = colorRampPalette(c("yellow", "darkblue"))(100))



