#===============================================================================
# Load libraries
#===============================================================================

library(Seurat)
library(ggplot2)
#install.packages('future')
library(future)
#===============================================================================
# Set working directory
#===============================================================================

setwd("~/CANCER RADIATION/Glioblastoma")

#===============================================================================
# Sample names and corresponding folders
#===============================================================================

sample.names <- c(
  "GBM131_Naive",
  "GBM131_RT2d",
  "GBM131_RT5d",
  "GBM131_RT3w",
  "GBM022_Naive",
  "GBM022_RT2d",
  "GBM827_Naive",
  "GBM827_RT2d"
)

sample.folders <- c(
  "GSM4967236",
  "GSM4967237",
  "GSM4967239",
  "GSM4967241",
  "GSM4967242",
  "GSM4967244",
  "GSM4967246",
  "GSM4967248"
)

#===============================================================================
# Read 10X data and create Seurat objects
#===============================================================================

seurat.list <- list()

for(i in 1:length(sample.names)){
  
  # Read 10X matrix
  counts <- Read10X(sample.folders[i])
  
  # Create Seurat object
  seurat.list[[sample.names[i]]] <- CreateSeuratObject(
    counts = counts,
    project = sample.names[i]
  )
  
}
#==============================================================================

#===============================================================================
# Calculate mitochondrial percentage
#===============================================================================
#Calculate mitochondrial percentage

for(i in 1:length(seurat.list)){
  
  seurat.list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    seurat.list[[i]],
    pattern = "^MT-"
  )
  
}

head(seurat.list[["GBM131_Naive"]]@meta.data)

names(seurat.list)

#===============================================================================
# QC: Violin plots
#===============================================================================

#===============================================================================
# Create QC directory
#===============================================================================

if(!dir.exists("QC")){
  dir.create("QC")
}

#===============================================================================
# QC: Violin plots and save as PDF
#===============================================================================

for(i in 1:length(seurat.list)){
  
  p <- VlnPlot(
    seurat.list[[i]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  ) +
    ggtitle(names(seurat.list)[i])
  
  print(p)
  
  ggsave(
    filename = paste0("QC/", names(seurat.list)[i], "_ViolinPlot.pdf"),
    plot = p,
    width = 10,
    height = 5
  )
  
}

#==============================================================================
# create scallter plot inside same qc directory
#===============================================================================

#===============================================================================
# QC: Scatter plots
#===============================================================================

if (!dir.exists('Scattered.plot')) {
  dir.create('Scattered.plot')
}


for(i in 1:length(seurat.list)){
  
  p1 <- FeatureScatter(
    seurat.list[[i]],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    
  )
  p2 <- FeatureScatter(
    seurat.list[[i]],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  p = p1 + p2
  
  print(p)
  
  ggsave(paste0("Scattered.plot/", names(seurat.list)[i],"scatteredPlot.pdf"),
         plot = p,
         height = 5,
         width = 10)
  
}

#==============================================================================
# filtering after QC
#==============================================================================

#===============================================================================
# Filter cells
#===============================================================================

for(i in 1:length(seurat.list)){
  
  seurat.list[[i]] <- subset(
    seurat.list[[i]],
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 6500 &
      nCount_RNA > 1000 &
      nCount_RNA < 45000 &
      percent.mt < 20
  )
  
}

#==============================================================================
# Calculate number of cells after filtering
#==============================================================================

for (i in 1:length(seurat.list)) {
  cat(names(seurat.list)[i], ":", ncol(seurat.list[[i]]),"Cells\n")
  
}

#===============================================================================
# lets save the filtered seurat list object before further downstream analysis
#==============================================================================


if (!dir.exists("filtered.object")) {
  dir.create("filtered.object")
  
}

for (i in 1:length(seurat.list)) {
  saveRDS(object = seurat.list[[i]],
    file =  paste0("filtered.object/",names(seurat.list)[i],"_RDS_filtered.rds"))
  
}
#===============================================================================
# change the future settings before calling SCTransform()
#===============================================================================


# Increase the maximum object size to 10 GB
options(future.globals.maxSize = 10 * 1024^3)

# Disable parallel processing
plan(sequential)
#===============================================================================
# SCTransform normalization
#===============================================================================

for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- SCTransform(seurat.list[[i]],
                                  vars.to.regress = "percent.mt",
                                  verbose = FALSE)
  gc()
  
}

#=============================================================================
#Phase 1 (Discovery):Use only the longitudinal series:

#===============================================================================
# Discovery cohort
#===============================================================================

gbm131.list <- seurat.list[c(
  "GBM131_Naive",
  "GBM131_RT2d",
  "GBM131_RT5d",
  "GBM131_RT3w"
)]

#===============================================================================
# Select integration features
#===============================================================================

features <- SelectIntegrationFeatures(
  object.list = gbm131.list,
  nfeatures = 3000
)

#===============================================================================
# Prepare SCT integration
#===============================================================================

gbm131.list <- PrepSCTIntegration(
  object.list = gbm131.list,
  anchor.features = features
)

#===============================================================================
# Find anchors
#===============================================================================

anchors <- FindIntegrationAnchors(
  object.list = gbm131.list,
  normalization.method = "SCT",
  anchor.features = features
)

#===============================================================================
# Merge
#===============================================================================

gbm131.merged <- merge(
  gbm131.list[[1]],
  y = gbm131.list[2:4],
  add.cell.ids = names(gbm131.list),
  project = "GBM131"
)
