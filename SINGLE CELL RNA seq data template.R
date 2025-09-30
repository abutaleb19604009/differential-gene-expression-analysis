# load the library
library(Seurat)
library(tidyverse)


# Load the NSCLC datasets

nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')

str(nsclc.sparse.m)

# extract the count matrix

cnt <- nsclc.sparse.m$`Gene Expression`
# build the seurat matrix
nsclc.seurat.obj <- CreateSeuratObject(
  counts = cnt,
  project = 'NSCLS',
  min.cells = 3,
  min.features = 200
)



# visualize the matrix
str(nsclc.seurat.obj)


# visualize the metadata and cheek for mitochondrial data

nsclc.seurat.obj@meta.data


# find the ^MT- sequence, mitochondrial sequence and add in a column of metadata


nsclc.seurat.obj[['percent.mt']] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = '^MT-')

# view metadata again

nsclc.seurat.obj@meta.data


# create violin plot to visualize scattered nature of metadata

VlnPlot(nsclc.seurat.obj, features = c('nCount_RNA','nFeature_RNA','percent.mt'))
FeatureScatter(nsclc.seurat.obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')+
  geom_smooth(method = 'lm')




# now filter the low quality cells based on nFeature_RNA and percent.mt



nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA>200 & nFeature_RNA<2500
                           & percent.mt<5)

# after filtering low quality cells, visualize the seurat object again
nsclc.seurat.obj

# normalization of data

nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# view seurat object after normalization


str(nsclc.seurat.obj)


# now identify the highly variable features or genes, because highly variable features are only 
#important for downstream analysis


nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = 'vst', nfeatures = 2000)


# identify top 10 highly variable features or genes 


top10 <- head(VariableFeatures(nsclc.seurat.obj),10)


# plot the variable feature
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)

LabelPoints(plot = plot1, points = top10, repel = TRUE)



# scalling data to reduce atrifects and unwanted variation 


all.gene <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.gene)


# perform linear dimentionality reduction use runpca to reduce herterogenicity by variable feature obj


nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize top 5 pca with top 5 genes 

print(nsclc.seurat.obj[['pca']], dims = 1:5, nfeatures = 5)



# visualizing the pca result with the heatmap 
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# determine the dimentionality of the data with elbowplot, and variance estimation 
ElbowPlot(nsclc.seurat.obj)

# after getting dimentionality, find neighbout and clustert the features

nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# find cluster with setted resolution


nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))

# view the metada for updated variable cluster
view(nsclc.seurat.obj@meta.data)

# view dimentionality plot with group by different column created in the metadata
# find the best resolution 

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)


# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)



# non-linear dimentionality reduction

reticulate::py_install(packages = 'umap-learn')

# run umap
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

# either set label = true, or Label cluster function to label individual cluster


DimPlot(nsclc.seurat.obj, reduction = 'umap')








