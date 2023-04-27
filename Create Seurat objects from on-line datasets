#################################################################################
# CREATE SEURAT OBJECT FROM EXPRESSION DATA AND METADATA 

## libraries
library(Matrix)
library(Seurat)

## Read Expression Data as full matrix
online_data = read.table("Data.csv", sep=",", header=T, row.names = 1)

## read metadata
## NOTE! Check separator from the metadata file! Otherwise error in `[.data.frame`(meta.data, , ii, drop = FALSE) (usually , or ;)
online_metadata = read.table("Metadata.csv", sep=",", header=T, row.names=1)

## Convert to sparse matrix
matrix_data = as(as.matrix(online_data), "sparseMatrix")

## Create Seurat object using sparse.mat and meta objects created above
seurat_object = CreateSeuratObject(matrix_data, assay = "RNA", min.cells = 5, min.features = 500, meta.data = online_metadata)


#################################################################################
# READ DATA FROM NCBI GEO (SOFT)

# source ("http://www.bioconductor.org/biocLite.R")
# https://www.nature.com/articles/s41467-018-06052-0.pdf?origin=ppub
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE...

library(Biobase)
library(GEOquery)
library(exprso)

# Open a GDS file
gse = getGEO(filename= 'family.soft.gz')

#view data
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]
names(GPLList(gse))
show(gse)
exprs(gse)

data = as.data.frame(exprs(gsm[[1]]))

# Read Expression Data as full matrix
online_data = read.table("Data.csv", sep=",", header=T, row.names = 1)

# read metadata
# NOTE! Check separator from the metadata file! Otherwise error in `[.data.frame`(meta.data, , ii, drop = FALSE) (usually , or ;)
online_metadata = read.table("Metadata.csv", sep=",", header=T, row.names=1)

# Convert to sparse matrix
matrix_data = as(as.matrix(online_data), "sparseMatrix")

# Create Seurat object using sparse.mat and meta objects created above
seurat_object = CreateSeuratObject(matrix_data, assay = "RNA", min.cells = 5, min.features = 500, meta.data = online_metadata)

