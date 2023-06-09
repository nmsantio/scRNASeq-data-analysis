#####################################################################################################################################
# LIBRARIES
library(Seurat)        # for seurat objects
library(dplyr)         # for seurat objects
library(data.table)    # for setDT function

# READ DATA 
EC_subset = readRDS(file = "Filepath\\EC_subset.rds")

#####################################################################################################################################

# IDENTIFY RELATIVE CANCER/CTRL FOLD CHANGES               
# determine the sample separator and column identities as above

# define the function to be used to detect the fold changes between cancer and control 
FoldChangeFunction = function(object, name) {
  identTest = paste(name, "Cancer", sep = "_");
  identControl = paste(name, "Ctrl", sep = "_");
  markers = FindMarkers(object, ident.1 = identTest, ident.2 =identControl, min.pct = FALSE, logfc.threshold = FALSE, verbose = TRUE);
  for (i in 1:length(colnames(markers))) {
    colnames(markers)[i] = paste(colnames(markers)[i], name, sep = "_")
  }
  return(markers);
}


# list the fold changes for all genes in each cluster (note: this will take some time!)
FoldChangeList = list(FoldChangeFunction(EC_subset, "Cluster1"),
                      FoldChangeFunction(EC_subset, "Cluster2"),
                      FoldChangeFunction(EC_subset, "Cluster3"),
                      FoldChangeFunction(EC_subset, "Cluster4"),
                      FoldChangeFunction(EC_subset, "Cluster5"),
                      FoldChangeFunction(EC_subset, "Cluster6"),
                      FoldChangeFunction(EC_subset, "Cluster7"),
                      FoldChangeFunction(EC_subset, "Cluster8"),
                      FoldChangeFunction(EC_subset, "Cluster9"))

# add gene names to first column with header gene_name and keep columns 1:3 (gene_name, p_val and avg_logFC) 
FoldChangeList_final = lapply(FoldChangeList, function(marker) {
  return(setDT(marker, keep.rownames = "gene_name")[,1:3]);
})

# combine the list of matrices (clusters) to one matrix according to gene names and save to working directory
FoldChangeMatrix_AllClusters = Reduce(function(...) merge(..., by="gene_name", all=TRUE), FoldChangeList_final)
write.table(FoldChangeMatrix_AllClusters, "Filename.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote= F)

# reorganize and rename: read data using gene names as row names
FoldChangeMatrix_AllClusters = read.table("Filename.txt", header = TRUE, sep = "", row.names = "gene_name")

# reorganize and rename: create a table with only fold changes (every other column, starting from column 2)
odd_indexes = seq(from = 2, to = length(FoldChangeMatrix_AllClusters), by = 2)
FoldChangesForHeatmaps = FoldChangeMatrix_AllClusters[, odd_indexes]

# reorganize and rename: simplify the column headers 
colnames(FoldChangesForHeatmaps)[1] = "Cluster1"
colnames(FoldChangesForHeatmaps)[2] = "Cluster2"
colnames(FoldChangesForHeatmaps)[3] = "Cluster3"
colnames(FoldChangesForHeatmaps)[4] = "Cluster4"
colnames(FoldChangesForHeatmaps)[5] = "Cluster5"
colnames(FoldChangesForHeatmaps)[6] = "Cluster6"
colnames(FoldChangesForHeatmaps)[7] = "Cluster7"
colnames(FoldChangesForHeatmaps)[8] = "Cluster8"
colnames(FoldChangesForHeatmaps)[9] = "Cluster9"

# save the fold change .txt file to working directory
write.table(FoldChangesForHeatmaps, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F)

#####################################################################################################################################

# PRINT RNA VALUE AVERAGES FOR BOTH SAMPLES                    

# determine the sample separator and column identities as above!
Split_seurat = SplitObject(EC_subset, split.by = "sampleID")
DefaultAssay(Split_seurat[["ctrl"]]) = "RNA" 
DefaultAssay(Split_seurat[["cancer"]]) = "RNA" 

# print average RNA gene expression per cluster/ sample, combine and save to working directory
ctrl_ave = AverageExpression(Split_seurat[["ctrl"]])
cancer_ave = AverageExpression(Split_seurat[["cancer"]])
ctrl_ave2 = ctrl_ave[["RNA"]]
cancer_ave2 = cancer_ave[["RNA"]]
ECave = cbind(ctrl_ave2, cancer_ave2)
write.table(ECave, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F) 

#####################################################################################################################################

# REFERENCES
# Satija lab (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)

                 
                   
                   
