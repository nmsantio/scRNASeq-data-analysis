####################################################################################################################################

# READ LIBRARIES
library(Seurat)
libray (dplyr)

####################################################################################################################################

# CREATE SEURAT OBJECT FROM 10X RAW DATA 
# (similar for each data set used for data integration)
data = Read10X(data.dir = "File path\\raw_feature_bc_matrix_filename")
seurat_object = CreateSeuratObject(counts = data, min.cells = 5, min.features = 500)    

# FILTER BASED ON RIBOSOMAL AND MITOCHONDRIAL GENE EXPRESSION
seurat_object = PercentageFeatureSet(seurat_object, pattern = "^mt-", col.name = "percent.mt")
seurat_object = PercentageFeatureSet(seurat_object, pattern = "^Rp[sl]", col.name = "percent.rb")
seurat_filtered = subset(seurat_object, subset = nFeature_RNA > 1500 & nFeature_RNA < 4500 & percent.mt < 3 & percent.rb < 15)

# NORMALIZE
seurat_normalized = NormalizeData(seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# FIND VARIABLE FEATURES
sample1 = FindVariableFeatures(seurat_normalized, selection.method = "vst", nfeatures = 5000)

####################################################################################################################################

# INTEGRATE DATA
sample1$sampleID = "ctrl"
sample2$sampleID = "cancer"
data_anchors = FindIntegrationAnchors(object.list = c(ctrl, cancer), anchor.features = 4500, dims = 1:30)
data_combined = IntegrateData(anchorset = data_anchors, dims = 1:30)
DefaultAssay(data_combined) = "integrated"
genes = rownames(data_combined)
seurat_scaled = ScaleData(data_combined, features = genes)   

####################################################################################################################################

# PRINCIPAL COMPONENT ANALYSIS (PCA) 
seurat_object_PCA = RunPCA(seurat_scaled, n_components = 0.95)
pca = seurat_object_PCA[["pca"]]
eigenValues = (pca@stdev)^2  
cumulative_PCA = cumsum(pca@stdev^2/sum(pca@stdev^2))
cumulative_PCA

# CLUSTER 
# dims and resolution depend on data
seurat_SNN = FindNeighbors(seurat_object_PCA, reduction = "pca", dims = 1:39)
seurat_cluster = FindClusters(seurat_SNN, resolution = 1) 

# TSNE  
seurat_TSNE = RunTSNE(seurat_cluster, reduction = "pca", dims = 1:39)
TSNEPlot(seurat_TSNE,  label = F, split.by = "sampleID", pt.size = 0.3)  

# UMAP              
seurat_UMAP = RunUMAP(seurat_cluster, reduction = "pca", dims = 1:39)
UMAPPlot(seurat_UMAP, label = F, split.by = "sampleID", pt.size = 0.3)

# SAVE
saveRDS(seurat_TSNE, file = "Filepath\\seurat_TSNE.rds")
saveRDS(seurat_UMAP, file = "Filepath\\seurat_UMAP.rds")

####################################################################################################################################

# REFERENCES
# Satija lab (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)

