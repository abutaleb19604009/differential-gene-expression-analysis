#Human Breast Cancer: Visium Fresh Frozen, Whole Transcriptome

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# load the data h5 in r

hdf5.obj<-  Read10X_h5(filename = 'Visium_Human_Breast_Cancer_raw_feature_bc_matrix.h5', 
                       use.names = TRUE, unique.features = TRUE)

# convert the data in seurat object

hbc.seurat <- CreateSeuratObject(hdf5.obj)

# view the metadata

hbc.seurat@meta.data

# there is no miochondrial info, all it

hbc.seurat[['percent.mt']] <- PercentageFeatureSet(hbc.seurat, pattern = '^MT-')


# visualize the metadata again

hbc.seurat@meta.data


#visualize scattered ness of metadata using violing plot

VlnPlot(hbc.seurat, features = c('nCount_RNA','nFeature_RNA','percent.mt'))



# visualize feature features relationship with scatter plot

p1<- FeatureScatter(hbc.seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')+
  geom_smooth(method = 'lm')

p2 <- FeatureScatter(hbc.seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mt')+
  geom_smooth(method = 'lm')

p1+p2


# normalization of data

hbc.seurat<- NormalizeData(hbc.seurat, normalization.method = 'LogNormalize', scale.factor = 10000)

# select highly varaible features and plot it

hbc.seurat <- FindVariableFeatures(hbc.seurat, selection.method = 'vst', nfeatures = 2000 )
# identify top highly variable gene 
top10 <- head(top10 <- head(VariableFeatures(hbc.seurat), 10))



# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hbc.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2


# scalling the data before performing linear dimentionality reduction

all.gene <- rownames(hbc.seurat)
hbc.seurat <- ScaleData(hbc.seurat, features = all.gene)

# perform linear dimentionality reduction to
hbc.seurat<- RunPCA(hbc.seurat, features = VariableFeatures(object = hbc.seurat))



# examine and visualize pca result with different ways


DimPlot(hbc.seurat, reduction = 'pca')

DimHeatmap(hbc.seurat, dims = 1:5, cells = 500)




# choosing true diemntionality of the dataset using elbowplot
ElbowPlot(hbc.seurat) # the plot shows that, top 15 pca are variable 



hbc.seurat <- FindNeighbors(hbc.seurat, dims = 1:15) #euclidean distance in PCA space
hbc.seurat <- FindClusters(hbc.seurat, resolution = 0.5)


# look at the cluster id of the first 10 cells
head(Idents(hbc.seurat),10)


# non- linear dimentionality reduction by UMAP AND tSNE

hbc.seurat <- RunUMAP(hbc.seurat, dims = 1:15)
DimPlot(hbc.seurat, reduction = 'umap')



# finding the differentially expressed features or cluster of biomarker 

# find all marker of cluster 2, all the gene differentially expressed in cluster2 cells

cluster2marker <- FindMarkers(hbc.seurat, ident.1 = 2)
head(cluster2marker, n=5)


# find all marker for every calster compared to remaining cells, report only +ve

everycluster.marker <- FindAllMarkers(hbc.seurat, only.pos = TRUE)

everycluster.marker %>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)


# chech the expression level of paricullar features across all the cluster
# distribution of feature across all cluster


VlnPlot(hbc.seurat, features = c('XBP1','FOXA1'))


# visualize with feature plot, expression of a feature across all cluster

FeaturePlot(hbc.seurat, features = c('XBP1','FOXA1'))

# plot heatmap with expression level across all cluster and include only top
# 20 features 
hbc.seurat.marker <- FindAllMarkers(hbc.seurat, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)

hbc.seurat.marker %>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC>1)%>%
  slice_head(n=10)%>%
  ungroup() ->top10
DoHeatmap(hbc.seurat, features = top10$gene)+ NoLegend()


FeaturePlot(hbc.seurat, features = c("EPCAM", "VIM", "PTPRC"))

# pathway enrichment 

install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))


# load the library
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


hbc.seurat.marker <- everycluster.marker

# Extract marker genes for cluster 2
cluster2.genes <- hbc.seurat.marker %>%
  dplyr::filter(cluster == 2 & avg_log2FC > 1) %>%
  pull(gene)


ego <- enrichGO(gene          = cluster2.genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",     # "BP" = Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)


# GO enrichment barplot
barplot(ego, showCategory = 10)

# KEGG enrichment dotplot
dotplot(ekegg, showCategory = 10)

# Enrichment map
emapplot(ego)

# Network plot of enriched pathways
cnetplot(ego, showCategory = 5)
