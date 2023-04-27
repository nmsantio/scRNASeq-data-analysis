####################################################################################################################################

# LIBRARIES
library (Seurat)
library (dplyr)
library(ggplot2)

####################################################################################################################################

# OPEN DATA

seurat_TSNE = readRDS(file = "Filepath\\seurat_TSNE.rds")
seurat_UMAP = readRDS(file = "Filepath\\seurat_UMAP.rds")

####################################################################################################################################

# EXAMPLE WITH THE UMAP FILE
Dataset = seurat_UMAP

# VISUALIZE GENE EXPRESSION IN DIFFERENT CLUSTERS TO DETERMINE CLUSTER IDENTITY

# PREPARE GENE LISTS AND PLOT BY DIFFERENT OPTIONS

# GENE LIST1                                 
Markers1 = c ("Pecam1", "Cdh5",                                           # endothelial
             "Rps12",                                                     # high ribo
             "Muc1",                                                      # epithelial
             "Pdgfa",                                                     # fibroblast
             "Pdgfrb", "Itga8",                                           # pericytes
             "Itga2b", "Itgb3",                                           # platelets
             "Ptprc",                                                     # white blood cells
             "Cd68",                                                      # macrophages                      
             "S100a8", "S100a9",                                          # neutrophils  
             "Ccr7",                                                      # lymphocytes
             "Cd28",                                                      # T cells
             "Bank1", "H2-Ob")                                            # B cells
                           
# SELECT DEFAULT ASSAY              
DefaultAssay(separatedCRECO) = "RNA"  

# DOTPLOT
DotPlot(Dataset, features = ((Markers1))) + RotatedAxis() + scale_colour_gradient2(low = "#0000FF", mid = "#FFFFFF", high = "#FF0000")

# VIOLIN PLOT
VlnPlot(Dataset, features = (Markers1), split.by = "sampleID", split.plot = TRUE) 

# FEATURE PLOT
FeaturePlot(Dataset, features = (Markers1), cols = c("grey", "red"), split.by = "sampleID")

# SAVE AS TIFF
tiff("Filename.tiff", units="in", width=8, height=3, res=300)
# code line here
dev.off()

# RENAME/SUBSET CLUSTERS
Subset1 = subset(Dataset, idents = c("Cluster1","Cluster2"))
Subset1_renamed = RenameIdents((Subset1),
                                 `Cluster1` = "Artery",
                                   `Cluster2` = "Vein")
                            
# SAVE/READ THE SEURAT OBJECT
saveRDS(Subset1_renamed, file = "Filepath\\Subset1_renamed.rds")
Subset1_renamed = readRDS(file = "Filepath\\Subset1_renamed.rds")

# PRINT CLUSTER CELL NUMBERS 
cluster_cellnro = table(Idents(Subset1_renamed), Subset1_renamed$sampleID)
write.table(cluster_cellnro, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F) 

# PRINT CLUSTER CELL PERCENTAGES
cluster_cellperc = prop.table(table(Idents(ClusterIdentities), Subset1_renamed$sampleID), margin = 2)
write.table(cluster_cellperc, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F) 

####################################################################################################################################

# RENAME
EC_Subset = Subset1_renamed

# IDENTIFY CLUSTER MARKERS 
ClusterMarkers = FindAllMarkers(EC_subset, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
write.table(ClusterMarkers, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F)

# IDENTIFY DIFFERENTIALLY EXPRESSED GENES IN ONE CLUSTER BETWEEN SAMPLES
# note: you can add more than two identities if you wish
# create a column in the meta.data slot to hold both cell type and sample ID information
# switch current identity to that column, find genes and save table to working directory
EC_subset$celltype_sampleID = paste(Idents(EC_subset), EC_subset$sampleID, sep = "_")
EC_subset$celltype = Idents(EC_subset)
Idents(EC_subset) = "celltype_sampleID"
DEgenes = FindMarkers(EC_subset, ident.1 = "Cancer", ident.2 ="Ctrl", min.pct = 0.25, logfc.threshold = 0.2)
write.table(DEgenes, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F)

####################################################################################################################################

# REFERENCES
# Satija lab (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)
