# reading single cell data into a seurat object in R

# the table is feature barcode sparce matrix
# the row is features or genes, the column is a cell

library(Seurat)
library(SeuratDisk)
library(loomR)
# read .RDS   FORMAT

rds.object <-readRDS('filename')

# check whether it is processed or raw rds object,bcz sometime processed data can also be read as rds



# 10x cell rangers .HDF5 FORMAT

hdf5.obj<-  Read10X_h5(filename = 'x', use.names = TRUE, unique.features = TRUE)

# this are just feature barcode matrix, we have to convert these into seurat object

seurat.hdf5 <- CreateSeuratObject(counts = hdf5.obj)



# .mtx file, if cell ranger analyze data, it produce mtx file, barcode file,and features file
# we can read all these into r and filnally convert into seurat object


mtx.obj <- ReadMtx(
  mtx = '',cells = '',features = ''
)

# then convert the feature barcode matrix into seurat object 

serut.mtx <- CreateSeuratObject(counts = mtx.obj)

# fix loom file 

library(loomR)
library(Seurat)

loom.obj <- connect("GSM4274191_CS20_longbone_rawdata.loom", mode = "r")

# Directly extract the expression matrix
expr_matrix <- loom.obj[["matrix"]][, ]

# Get gene and cell names
genes <- loom.obj$row_attrs$Gene[]
cells <- loom.obj$col_attrs$CellID[]

rownames(expr_matrix) <- genes
colnames(expr_matrix) <- cells

# Create Seurat object manually
loom2seurat <- CreateSeuratObject(counts = expr_matrix)

loom.obj$close_all()

loom2seurat



# read the .loom file 

loom.obj <- connect(filename = 'GSM4274191_CS20_longbone_rawdata.loom', mode = 'r')

loom.obj
# convert the loom object into a seurat object

loom2seurat <- as.Seurat(loom.obj)



