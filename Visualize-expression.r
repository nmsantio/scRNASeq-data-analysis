########################################################################################################################
# LIBRARIES
library(Seurat)        
library(heatmaply)    
library(dplyr)        
library(RColorBrewer)  
library(stringr)       
library(plotly)
library(htmltools)

#######################################################################################################################
# READ FOLD CHANGES OR RNA DATA FROM WORKING DIRECTORY 

# fold changes (FCs) or average RNA expression (RNAave)
# change column order if needed and save
FCs = read.table("Filename.txt", sep="\t", header = TRUE)
Ordered = data.frame(FCs[,c(6,1,5,3,4,2,7)])
write.table(Ordered, "FCs_filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F) 
FCs = read.table("FCs_filename.txt", sep="\t", header = TRUE)

RNAave = read.table("Filename.txt", sep="\t", header = TRUE)
Ordered = data.frame(RNAave[,c(6,1,5,3,4,2,7)])
write.table(Ordered, "RNAave_filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F) 
RNAave = read.table("RNAave_filename.txt", sep="\t", header = TRUE)

#####################################################################################################################
# ADD GENES MANUALLY 

genelist1 = c("Vegfa", "Vegfb", "Vegfc", "Vegfd", "Flt1", "Kdr", "Flt4")

# GENE LIST FROM DATABASE
# save as txt and open, select gene names and rename the column 
genelist2 = data.frame (Filename$Column)
# remove duplicates and rename
genelist2 = unique(genelist)
colnames(genelist) = c("Gene")
# add quotes ("")
genelist2 = toString(shQuote(Filename$Column))
# save to working directory
write.table(genelist2, "genelist2.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote= F)

#####################################################################################################################
#  PREPARE HEATMAPS FROM FOLD CHANGES AND RNA DATA 

# FC data to heatmap without dendrogram 
hitset = FCs[c(genelist1),]                                    # subset hits from cancer/ctrl FC txt file
expressionmatrix = data.matrix ((hitset))                           # data.frame to matrix 
colBWR = colorRampPalette(c("blue", "white", "red"))                # color palette BWR
heatmaply(expressionmatrix, dendrogram = "none", colors = colBWR,  
          limits = c(-1, 1))

# expression data to heatmap 
hitset = RNAave[genelist1,]                                           # subset hits from RNA average txt file
expressionmatrix = data.matrix ((hitset))                           # data.frame to matrix 
BuYlRd = colorRampPalette(c("#4475B4","#73ACD1","#ABD9E9",          # manually listing RdYlBu colors in reverse (BuYlRd)
                            "#E0F3F8","#FFFDBF","#FEE090",
                            "#FDAE61","#F46D43", "#D73027" ))    
heatmaply(expressionmatrix, dendrogram = "none", colors = BuYlRd)   # heatmap without dendrogram

#extract values from the heatmap and save as txt
heatmapvalues = subset(FCs, rownames(FCs) %in% genelist1)
write.table(heatmapvalues, "Filename.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote= F)

#####################################################################################################################
# VISUALIZE THE CLUSTERS BY FEATURE PLOTS             

# read tSNE or UMAP Seurat object
# set default and identify sample IDs
EC_subset = readRDS(file = "Filepath\\filename.rds")a
DefaultAssay(EC_subset) = "RNA"       
objects = SplitObject(EC_subset, split.by = "sampleID")
DefaultAssay(objects[["ctrl"]]) = "RNA" 
DefaultAssay(objects[["cancer"]]) = "RNA" 

# visualize by feature plots (genes, features etc. replace names as you wish)
FeaturePlot(objects[["ctrl"]], features = c("Gene1", "Gene2", "Gene3")) 
FeaturePlot(objects[["cancer"]], features = c("Gene1", "Gene2", "Gene3")) 

# save to working directory (replace the plot code depending on what you wish to save)
tiff("Filename.tiff", units="in", width=8, height=7, res=300)
FeaturePlot(objects[["ctrl"]], features = c("Gene1", "Gene2", "Gene3"), max.cutoff = 3, pt.size = 0.4, cols = c("grey", "red"))
dev.off()

# plot and save as svg for figure preparation
svg("Filename.svg")
FeaturePlot(objects[["ctrl"]], features = c("Gene1", "Gene2", "Gene3"), pt.size = 0.2, max.cutoff = 3)
dev.off()

####################################################################################################################
# VISUALIZE INTEGRATED DATA BY DOHEATMAP

DefaultAssay(EC_subset) = "integrated"   
DoHeatmap(EC_subset, features = rev(genelist1)) +  scale_fill_gradientn(colors = c("blue", "white", "red") )

####################################################################################################################
# REFERENCES
# Satija lab (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)

         


