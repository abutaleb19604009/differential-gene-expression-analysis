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

rds.object@meta.data$nCount_RNA
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

rds.object <- RunUMAP(rds.object, dims = 1:25)
rds.object <- RunTSNE(rds.object, check_duplicates = FALSE )

# find neighbout and cluster

rds.object <- FindNeighbors(rds.object, dims = 1:15)

rds.object <- FindClusters(rds.object, resolution = 0.3)

DimPlot(rds.object, reduction = 'umap')


# Single r cell type annotation 

# convert the seurat object into single cell experiment 
BiocManager::install("SingleCellExperiment")

# as celldex is a part of Rhdf5lib, we need to first install this then celldex
BiocManager::install("Rhdf5lib")
BiocManager::install("celldex")

BiocManager::install("SingleR")
library(SingleR)
library(celldex)
library(SingleCellExperiment)

SCE <- Seurat::as.SingleCellExperiment(rds.object)

# load the reference datasets

ref <- celldex::HumanPrimaryCellAtlasData()


# run single cell annotation

pred_HCA <- SingleR::SingleR(test = SCE, ref = ref, assay.type.test=1,
                             labels = ref$label.main)


# convert to dataframe 

pred_HCA <- as.data.frame(pred_HCA)
pred_HCA$Barcode <- rownames(pred_HCA)

# now mapped to seurat cells

rds.object$SingleR_label <- pred_HCA$labels[
  match(colnames(rds.object), pred_HCA$Barcode)
]




# now visualize the result

library(Seurat)

DimPlot(rds.object, group.by = "SingleR_label", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotation by SingleR (HumanPrimaryCellAtlas)")
