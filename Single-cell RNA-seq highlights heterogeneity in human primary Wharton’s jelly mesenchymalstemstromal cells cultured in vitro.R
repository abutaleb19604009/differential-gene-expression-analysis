# reading single cell data into a seurat object in R

# the table is feature barcode sparce matrix
# the row is features or genes, the column is a cell

library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)



# read .RDS   FORMAT

rds.object <-readRDS('1GSE220994_AF_scRNAseq_seurat_object.rds')

rds.object@meta.data
rds.object@assays$RNA
rds.object[['RNA']]
rds.object@reductions$pca %>% head()
rds.object[['pca']]




#visualize scattered ness of metadata using violing plot

VlnPlot(rds.object, features = c('nCount_RNA','nFeature_RNA','percent.mt'))

# check normalization was done or not
DefaultAssay(rds.object) <-"RNA" 

rds.object[["RNA"]]@data[1:5,1:5]



# find variable features was found or not
length(VariableFeatures(rds.object)) # if it is 2000, then variable features was found, if it is zero, then it was not found

head(VariableFeatures(rds.object))

# check scalling was done or not. if it is 0x0 matrix then scalling is not done

rds.object[["RNA"]]@scale.data


#check about pca, if it is showing null, pca is not performed

rds.object@reductions$pca



# as data is not normalized, this pca is wrong, remove current pca
rds.object@reductions$pca <- NULL

# set the default assayto rna

DefaultAssay(rds.object) <- "RNA"

# normalize data

rds.object <- NormalizeData(rds.object)


# find variable features


rds.object <- FindVariableFeatures(rds.object)


# scale the data

rds.object <- ScaleData(rds.object)

# runpca
rds.object <- RunPCA(rds.object)

DimHeatmap(rds.object, dims = 1:5, reduction = 'pca', cells = 500)

#visualizing pca
ElbowPlot(rds.object) # 15 pca has choosen 

# run umap

rds.object <- RunUMAP(rds.object, dims = 1:15)
rds.object <- RunTSNE(rds.object, check_duplicates = FALSE )

# find neighbout and cluster

rds.object <- FindNeighbors(rds.object, dims = 1:15)

rds.object <- FindClusters(rds.object, resolution = 0.3)

DimPlot(rds.object, reduction = 'umap')

# load the metadata csv file

meta <- read.csv("GSE220994_AF_scRNAseq_metadata.csv", stringsAsFactors = FALSE)

# set the rowname to barcode

rownames(meta) <- meta$Barcode

#metadata rows must be in the same order as colnames(rds.object)

meta <- meta[colnames(rds.object), ]

# check again they are smae order or not

head(colnames(rds.object))
head(rownames(meta))

# add meta data to seurat object
rds.object <- AddMetaData(rds.object, metadata = meta, )

# now check the colname of rds object metadata and see trimester information is 
# there or not 

colnames(rds.object@meta.data) # trimester is added
table(rds.object$Trimester)

# now dimplot visualize cells by trimester 

DimPlot(rds.object, reduction = 'umap', group.by = 'Trimester')

# subset seurat object based on trimester 

T2 <- subset(rds.object, subset = Trimester == "Second")
T3 <- subset(rds.object, subset = Trimester == "Third")

#split by Trimester in plots
DimPlot(rds.object, reduction = "umap", split.by = "Trimester")
DimPlot(rds.object, reduction = "umap", label = TRUE)

# now annotate the cluster manually, for stem cell automatic annotation is not 
# possible 

marker <- FindAllMarkers(
  rds.object,
  min.pct = 0.25,
  only.pos = TRUE,
  logfc.threshold = 0.25
)

# annotate the cells
BiocManager::install("singleR")
library(SingleR)
library(celldex)

# convert seurat object into SingleCellExperiment 
library(SingleCellExperiment)

single_cell_experiment <- as.SingleCellExperiment(rds.object)

# load human reference cells data sets
ref <- celldex::HumanPrimaryCellAtlasData()

# run automatic cell type annotation with singleR

pred <- SingleR(
  test = single_cell_experiment,
  ref = ref,
  labels = ref$label.main
)

# add annotated cells into seurat object 

rds.object$CellType_SingleR <- pred$labels

# visualize annotated cells using groupby 

DimPlot(rds.object, group.by = "CellType_SingleR" )
DimPlot(rds.object, reduction = 'umap', group.by = 'CellType_SingleR', 
        split.by ="Trimester" )
