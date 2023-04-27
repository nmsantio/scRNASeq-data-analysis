# LIBRARIES
library(slingshot)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(broom)
library(Seurat)

# UMAP SEURAT OBJECT
# create a column in the meta.data slot to hold cell type information
UMAP = readRDS(file = "Filepath\\Subset1_renamed.rds")
UMAP$celltype = paste(Idents(UMAP))

# TRAJECTORY ANALYSIS BY SLINGSHOT
sds = slingshot(Embeddings(UMAP, "umap"), clusterLabels = UMAP$celltype, 
                 stretch = 0)

cell_pal = function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal = pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal = setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors = cell_pal(UMAP$celltype, brewer_pal("qual", "Set4"))
cell_colors_clust = cell_pal(UMAP$seurat_clusters, hue_pal())

tiff("Filename.tiff", units="in", width=6, height=6, res=300)
plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.3)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()

# REFERENCE TO SLINGSHOT
# https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html 


