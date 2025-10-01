# load the library
library(tidyverse)
library(ggplot2)
library(Seurat)
library(gridExtra)

# list the directory

dirs <- list.dirs(path ='ssc.h.c/', recursive = F, full.names = F )

# 

for (x in dirs) {
  name<- gsub('_filtered_feature_bc_matrix','',x)
  cnt<- ReadMtx( mtx = paste0('data/',x,'/matrix.mtx.gz'),
           features = paste0('data/',x,'/features.tsv.gz'),
           cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cnt))
}

# run ls() for finding out oll the object in directory 
# merge the data into one seurat object so that quality scoring can be performed in once


metged.seurat <- merge(HB17_background, y= c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,HB53_background,
                            HB53_tumor), add.cell.ids= ls()[3:9],project= 'HB')
      
# view the metadata of metged.seurat
metged.seurat@meta.data

# create column in metadata
metged.seurat@meta.data # see the metadta and decide to copy row name of metadata as a colum

metged.seurat$sample <- rownames(metged.seurat@meta.data)

# split the new column into three more column



metged.seurat@meta.data<- separate(metged.seurat@meta.data, col = 'sample', into = c('patient','type','barcode'), sep = '_')

#view metadata again 

metged.seurat@meta.data




# calculate mitochondrial gene percentage

metged.seurat$percent.mt <- PercentageFeatureSet(metged.seurat, pattern = '^MT-')

#view metadata again


metged.seurat@meta.data


# applying QC on the combined Seurat object

# step1: filtering the cell 


metged.seurat.filtered <- subset(metged.seurat, subset = nCount_RNA >800 & nFeature_RNA > 500 &
                                   percent.mt <10)
# again view filtered data after filtering, 
metged.seurat.filtered #67851 cell are there

metged.seurat # before filtering it was 77936 


# perform standard workflow to find out if there is batch effects or not

# step2 normalization 

metged.seurat.filtered <- NormalizeData(object =  metged.seurat.filtered)

# step3 find variable features, only these are meaningful in downstream analysis

metged.seurat.filtered <- FindVariableFeatures(object= metged.seurat.filtered)

# step4 scale the data
metged.seurat.filtered <- ScaleData(object = metged.seurat.filtered)

# step 5 perform linear dimentionality reduction after scalling 

metged.seurat.filtered <- RunPCA(object = metged.seurat.filtered)

#step 6 elboplot, find the dimention of dataset after linear dimentionality reduction 

ElbowPlot(metged.seurat.filtered)
# 1:20 principle componenet has variation

#step7 find neighbour

metged.seurat.filtered <- FindNeighbors(object = metged.seurat.filtered, dims = 1:20)

# step8 find cluster 

metged.seurat.filtered <- FindClusters(object = metged.seurat.filtered)

# step run umap

metged.seurat.filtered <- RunUMAP(object = metged.seurat.filtered, dims = 1:20)

# now visualize the results with the dimplot

p1<- DimPlot(metged.seurat.filtered, reduction = 'umap', group.by = 'patient')
p2<- DimPlot(metged.seurat.filtered, reduction = 'umap', group.by = 'type', cols = c('red','blue','green'))

grid.arrange(p1,p2,nrow =2, ncol=2)

# as batch effects is coming due to the grouping data by patient
# split object should perform along patient column to reduce batch effect


obj.list <- SplitObject(metged.seurat.filtered, split.by = 'patient')

# for eacch obj in the list (each patient), do normalization and identify highly variable features


for (i in 1:length(obj.list) ) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  
}

# select integration features from the obj.list

featuress <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchore using selected features

anchores <- FindIntegrationAnchors(object.list = obj.list, anchor.features = featuress)

# build the integrated data

seurat.integrated <- IntegrateData(anchorset = anchores) 

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
