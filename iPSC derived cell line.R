

# load the library 
library(Seurat)

# see the all files present in workig directory and see if they are still in
# .gz format or not if yes then open linux and do gunzip 
list.files("C:/Users/pc/Documents/RNA_seq/mesenchymal stromal cells and stem cells by single-cell transcriptomic analysis/iPSC lines from the American Type Culture/")


# create count matrix for all the sample 
countDI <- ReadMtx(
  mtx = 'GSE208625_D_I_filtered_matrix.mtx',
  features = 'GSE208625_D_I_filtered_features.tsv',
  cells = 'GSE208625_D_I_filtered_barcodes.tsv'
)

countDN <- ReadMtx(
  mtx = 'GSE208625_D_N_filtered_matrix.mtx',
  features = 'GSE208625_D_N_filtered_features.tsv',
  cells = 'GSE208625_D_N_filtered_barcodes.tsv'
)

countNI <- ReadMtx(
  mtx = 'GSE208625_N_I_filtered_matrix.mtx',
  features = 'GSE208625_N_I_filtered_features.tsv',
  cells = 'GSE208625_N_I_filtered_barcodes.tsv'
)

countNN <- ReadMtx(
  mtx = 'GSE208625_N_N_filtered_matrix.mtx',
  features = 'GSE208625_N_N_filtered_features.tsv',
  cells = 'GSE208625_N_N_filtered_barcodes.tsv'
)

# create seurat object for all these count matrix

DI <- CreateSeuratObject(counts = countDI)
DN <- CreateSeuratObject(counts = countDN)
NI <- CreateSeuratObject(counts = countNI)
NN <- CreateSeuratObject(counts = countNN)

# insert sample coulumn in metadata to seperate 
DI$sample <- 'DI'
DN$sample <- 'DN'
NI$sample <- 'NI'
NN$sample <- 'NN'

# views metadata 


DI@meta.data

# make a combined seurat object now
combined.seurat <- merge(DI, y=list(DN,NI,NN), 
                         add.cell.ids =c('DI','DN','NI','NN'))

# follow standard workflow 

combined.seurat$mt <- PercentageFeatureSet(combined.seurat, pattern = "^MT")
combined.seurat@meta.data

# visualize the parameter of metadata and identify cells with very high 
# and very low features, counts and mt should be less than 20
VlnPlot(combined.seurat, features = c('nFeature_RNA','nCount_RNA','mt'),
        ncol = 3)


combined.filtered <- subset(
  combined.seurat,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 7500 &
    nCount_RNA > 1000 &
    nCount_RNA < 100000 &
    mt < 20
)


# then follow standard workflow NFS PCA NEIG CLUSTER UMAP

combined.filtered <- NormalizeData(combined.filtered)

combined.filtered <- FindVariableFeatures(combined.filtered)

combined.filtered <- ScaleData(combined.filtered)


combined.filtered <- RunPCA(combined.filtered)
ElbowPlot(combined.filtered) # dim = 1:20


combined.filtered <- FindNeighbors(combined.filtered)
combined.filtered <- FindClusters(combined.filtered)


combined.filtered <- RunUMAP(combined.filtered, dims = 1:20)




DimPlot(combined.filtered, reduction = 'pca')

DimPlot(combined.filtered, reduction = 'umap')

DimPlot(combined.filtered, reduction = 'umap', group.by = 'sample')


# the dimplot shows that cluster is due to batch effects of sample, not 
# due to real biological differences 


list.objt <- SplitObject(combined.filtered, split.by = 'sample')

for (i in 1:length(list.objt)){
  list.objt[[i]] <- NormalizeData(list.objt[[i]])
  list.objt[[i]] <- FindVariableFeatures(list.objt[[i]])
}


# now find integration features and integration anchores using selected features

featuress <- SelectIntegrationFeatures(object.list = list.objt)
anchorage <- FindIntegrationAnchors(object.list = list.objt, anchor.features = 
                                      featuress)


# now build the integrated data

integrated.seurat <- IntegrateData(anchorset = anchorage)


# scale runpca runumap and dimplot

integrated.seurat <- ScaleData(integrated.seurat)

integrated.seurat <- RunPCA(integrated.seurat)

integrated.seurat <- RunUMAP(integrated.seurat)



# now visualize againg with umap to see batch effect still persist or not 



DimPlot(integrated.seurat, group.by = 'sample')

DimPlot(integrated.seurat, group.by = '')
