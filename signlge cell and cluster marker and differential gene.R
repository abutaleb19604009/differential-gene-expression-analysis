library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# create a seurat object 

pbmc.seurat <- CreateSeuratObject(counts = pbmc.data, project = 'PBMC', min.cells = 3, min.features = 300)

# view metadata

pbmc.seurat@meta.data

# calculate percent.mt sequence

pbmc.seurat[['percent.mt']] <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')


# revisualize the meatadata

pbmc.seurat@meta.data



#visualize scattered ness of metadata using violing plot

VlnPlot(pbmc.seurat, features = c('nCount_RNA','nFeature_RNA','percent.mt'))

# visualize feature features relationship with scatter plot

p1<- FeatureScatter(pbmc.seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')+
  geom_smooth(method = 'lm')

p2 <- FeatureScatter(pbmc.seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mt')+
  geom_smooth(method = 'lm')

p1+p2


# normalization of data

pbmc.seurat<- NormalizeData(pbmc.seurat, normalization.method = 'LogNormalize', scale.factor = 10000)

# select highly varaible features and plot it

pbmc.seurat <- FindVariableFeatures(pbmc.seurat, selection.method = 'vst', nfeatures = 2000 )
# identify top highly variable gene 
top10 <- head(top10 <- head(VariableFeatures(pbmc.seurat), 10))

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1  plot2


# scalling the data before performing linear dimentionality reduction

all.gene <- rownames(pbmc.seurat)
pbmc.seurat <- ScaleData(pbmc.seurat, features = all.gene)


# perform linear dimentionality reduction to
pbmc.seurat<- RunPCA(pbmc.seurat, features = VariableFeatures(object = pbmc.seurat))

# examine and visualize pca result with different ways


DimPlot(pbmc.seurat, reduction = 'pca')

DimHeatmap(pbmc.seurat, dims = 1:5, cells = 500)


# choosing true diemntionality of the dataset using elbowplot
ElbowPlot(pbmc.seurat) # the plot shows that, top 12 pca are variable 



# find neighbourhood and find cluster

pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:12) #euclidean distance in PCA space
pbmc.seurat <- FindClusters(pbmc.seurat, resolution = 0.5)

# look at the cluster id of the first 10 cells
head(Idents(pbmc.seurat),10)

# non- linear dimentionality reduction by UMAP AND tSNE

pbmc.seurat <- RunUMAP(pbmc.seurat, dims = 1:12)
DimPlot(pbmc.seurat, reduction = 'umap')


# finding the differentially expressed features or cluster of biomarker 

# find all marker of cluster 2, all the gene differentially expressed in cluster2 cells

cluster2marker <- FindMarkers(pbmc.seurat, ident.1 = 2)
head(cluster2marker, n=5)


# find all the gene distinguish cluster 2 and cluster 3

cluster2.3marker <- FindMarkers(pbmc.seurat, ident.1 = 2, ident.2 = 3)

head(cluster2.3marker, n=5)


# find all the gene distinguish cluster 2 and cluster 3 &4

cluster2.34marker <- FindMarkers(pbmc.seurat, ident.1 = 2, ident.2 = c(3,4))

head(cluster2.3marker, n=5)

# find all marker for every calster compared to remaining cells, report only +ve

everycluster.marker <- FindAllMarkers(pbmc.seurat, only.pos = TRUE)

everycluster.marker %>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)



# find claster0 marker against all other cells
# it select the cells with logfc >=0.25 and with statistical test =
#ROC test returns the ‘classification power’ for any individual marker 
#(ranging from 0 - random, to 1 - perfect).

claster0marker <- FindMarkers(pbmc.seurat, ident.1 = 0,logfc.threshold = 0.25, 
                              test.use = "roc", only.pos = TRUE)

claster0marker



# chech the expression level of paricullar features across all the cluster
# distribution of feature across all cluster


VlnPlot(pbmc.seurat, features = c('RPS6','RPS25'))


# we can also plot raw count in vlnplot

VlnPlot(pbmc.seurat, features =c('RPS6','RPS25'), slot = 'counts', log = TRUE )



# visualize with feature plot, expression of a feature across all cluster

FeaturePlot(pbmc.seurat, features = c('RPS6','RPS25','PEBP1'))




# plot heatmap with expression level across all cluster and include only top
# 20 features 
pbmc.seurat.marker <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)

pbmc.seurat.marker %>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)%>%
  slice_head(n=10)%>%
  ungroup() ->top10
DoHeatmap(pbmc.seurat, features = top10$gene)+ NoLegend()

# first find all marker for each cluster
pbmc.markers <- FindAllMarkers(pbmc.seurat, only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
#search databases (like PanglaoDB or CellMarker)to convert marker into cells 
# identity


# find out level of cluster or based on which gene cluster is created
# is nothing but level of seurat object


levels(pbmc.seurat)
# assigning cell type identity to cluster

new.cell.identity <- c( c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                      "NK", "DC", "Platelet"))

names(new.cell.identity) <- levels(pbmc.seurat)

pbmc.seurat<- RenameIdents(pbmc.seurat, new.cell.identity)

DimPlot(pbmc.seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
