# Data: GSE157657, Last Modification date: 11/03/2026
#==============================
# Load libraries
# ===============================
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
library(limma)
install.packages("pheatmap")
library(pheatmap)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(dplyr)

setwd("C:/Users/abuta/Downloads")
# ===============================
# Load expression data
# ===============================
counts <- read.table(
  "norm.data.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# remove the gene nane and only keep the expression matrix
expr <- as.matrix(counts[, -1])

# keep only gene name column for gene annotation
gene_annot <- counts[, 1, drop = FALSE]

# ===============================
# Load metadata
# ===============================

gse <- getGEO("GSE157657", GSEMatrix = TRUE)
meta <- pData(gse[[1]])
write.csv(meta, "metadata.csv")

# KEY FIX: use internal sample IDs
rownames(meta) <- meta$title

# ===============================
# Align metadata to expression, see col name of expression matrix match row name of meta
# keep new metadata
# ===============================
meta <- meta[colnames(expr), ]

# FINAL sanity check (must be TRUE)
stopifnot(all(colnames(expr) == rownames(meta)))


all(colnames(expr) == rownames(meta))

#======================================

# take only necessary column from the metadata and clean the column name and content


meta.clean <- meta %>%
  select(title,
         geo_accession,
         characteristics_ch1,
         characteristics_ch1.1,
         characteristics_ch1.2,
         characteristics_ch1.3,
         characteristics_ch1.5,
         characteristics_ch1.6,
         characteristics_ch1.7)

#rename the column name of the meta.clean metadata

meta.clean <- meta.clean %>%
  rename(
    patient_id = characteristics_ch1,
    subgroup = characteristics_ch1.1,
    group = characteristics_ch1.2,
    days_from_att = characteristics_ch1.3,
    subgroup_att = characteristics_ch1.5,
    response_group = characteristics_ch1.6,
    gender = characteristics_ch1.7
  )

# now mutate the value of each variable 


meta.clean <- meta.clean %>%
  mutate(
    patient_id = gsub("patient id:","",patient_id),
    subgroup = gsub("subgroup:","",subgroup),
    group = gsub("group:","",group),
    days_from_att = as.numeric(gsub("days_from_att:","",days_from_att)),
    subgroup_att = gsub("subgroup_att:","",subgroup_att),
    response_group = gsub("response_group:","",response_group),
    gender = gsub(":","", gender)
    
    
    
  )
write.csv(meta.clean,"meta.clean.inconsistant.csv" ,row.names = TRUE) 
meta.consistant.clean<- read.csv("meta.consistant .csv")

#======================================================================


# force subgroup to clean character
meta.clean$subgroup <- as.character(meta.clean$subgroup)
meta.clean$subgroup <- trimws(meta.clean$subgroup)
meta.clean$subgroup <- gsub("\\s+", " ", meta.clean$subgroup)

# now refactor
meta.clean$subgroup <- factor(meta.clean$subgroup)
  

#=====================================================================

# analyse the data based on subgroup column of meta
unique(meta.clean$subgroup)

# check sample distribution of each group 
table(meta.clean$subgroup)

#============================================================================

# create a subset of meta.clean and expr for limma analysis by subgroup
# 1st analysis for drug resistance vs control

keep1_DR_C <- meta.clean$subgroup %in% c("TB Drug Resistance","Control")


meta.clean.DR_C <- meta.clean[keep1_DR_C,]
expr.sub1_DR_C <- expr[,keep1_DR_C]

# now convert the subgroup in the meta.clean.sub1 into factor

meta.clean.DR_C$subgroup <- factor(
  meta.clean.DR_C$subgroup,
  levels = c("Control", "TB Drug Resistance"))
  all(colnames(expr.sub1_DR_C) == rownames(meta.clean.DR_C))


#===========================================================================

  # now run limma for the analysis of "Control" vs "TB Drug Resistance"
  library(limma)
  
  # Make sure subgroup is a factor with correct reference
  
  # 
  # meta.clean.sub1$subgroup <- factor(
  #   meta.clean.sub1$subgroup,
  #   levels = c("Control", "TB Drug Resistance")
  # )
  
  # Design matrix
  design <- model.matrix(~ subgroup, data = meta.clean.DR_C)
  
  # Fit model
  fit <- lmFit(expr.sub1_DR_C, design)
  fit <- eBayes(fit)
  
  # Get results
  res <- topTable(
    fit,
    coef = "subgroupTB Drug Resistance",
    number = Inf,
    adjust.method = "BH"
  )
  
  head(res)
  
table(head(res))


#==============================================================================

# now filter out differentially expressed genes
library(tibble)
deg <- res %>%
  rownames_to_column("Ensembl_ID") %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(deg)

#===============================================================

# represent data by preapring a volcano plot

library(EnhancedVolcano)

p<-EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "logFC",
  y = "adj.P.Val",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TB Drug Resistance vs Control",
  
)

#=====================================================
#save the generated plot

ggsave(
  filename = "TB_DrugResistance_vs_Control_volcano.png",
  plot = p,
  width = 8,
  height = 8,
  dpi = 300
)

#==========================================================================

# generate final DEG table and save it

deg <- res %>%
  rownames_to_column("Ensembl_ID") %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1)

write.csv(deg, "DEG_TB_DrugResistance_vs_Control.csv", row.names = FALSE)



# spli upregulated and downregulated genes

deg_up   <- deg %>% filter(logFC > 0)
deg_down <- deg %>% filter(logFC < 0)


#============================================================================

#pathway and functional enrichment 

BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

gene_ids <- deg$Ensembl_ID

ego <- enrichGO(
  gene          = gene_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dplot<- dotplot(ego, showCategory = 20)

ggsave(filename = "dot_plot_control_vs_drugresistant.jpg",
       plot = dplot,
       height = 8,
       width = 8,
       dpi = 300)

#=============================================================================
# generate heatmap for better graphical representation

# build the annotation column for the generation of heatmap
annotation_col <- data.frame(
  subgroup = meta.clean.sub1$subgroup
)

# THIS LINE IS CRITICAL
rownames(annotation_col) <- rownames(meta.clean.sub1)

# sanity check
all(rownames(annotation_col) == colnames(expr.sub1))



# now gererate pheatmap

library(pheatmap)

pheatmap<-pheatmap(
  expr.sub1[top_genes, ],
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE
)

ggsave(filename = "pheatmap_control_vs_drugresistant.jpg",
       plot =pheatmap,
       width = 8,
       height = 8,
       dpi = 300)

#========================================================================

# find out upregulated or downregualted pathway
# for upregulated pathway


ego_up <- enrichGO(
  gene          = deg_up$Ensembl_ID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot_up <- dotplot(ego_up, showCategory = 15, title = "Upregulated Pathways")

ggsave(
  "DR_Cimmune_upregulated_pathways.jpg",
  plot = dotplot_up,
  width = 8,
  height = 8,
  dpi = 300
)

# now for down regulation 


ego_down <- enrichGO(
  gene          = deg_down$Ensembl_ID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)


dotplot_down <- dotplot(ego_down, showCategory = 15, title = "Downregulated Pathways")

ggsave(
  "DR_Cimmune_downregulated_pathways.jpg",
  plot = dotplot_down,
  width = 8,
  height = 8,
  dpi = 300
)

# focus only immune pathways

immune_up <- ego_up@result %>%
  dplyr::filter(grepl("immune|interferon|cytokine|inflammatory|T cell|B cell|macrophage|neutrophil",
                      Description,
                      ignore.case = TRUE))

head(immune_up)


write.csv(
  immune_up,
  "immune_pathways_upregulated_TB_DR_vs_Control.csv",
  row.names = FALSE
)


# create dot plot only for immune pathway


library(ggplot2)

immune_plot <- ggplot(
  immune_up[1:15, ],
  aes(x = GeneRatio, y = reorder(Description, GeneRatio))
) +
  geom_point(aes(size = Count, color = p.adjust)) +
  labs(
    title = "Upregulated Immune Pathways",
    x = "Gene Ratio",
    y = "Immune Pathway"
  ) +
  theme_bw()

immune_plot

# save dot plot only for immune pathway


ggsave(
  "DR_Cimmune_signature_upregulated.jpg",
  plot = immune_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Enrichment plot

library(enrichplot)

ego_up <- pairwise_termsim(ego_up)

emap <- emapplot(ego_up)

emap
ggsave(
  "DR_Cimmune_enrichment_network.jpg",
  plot = emap,
  width = 8,
  height = 6,
  dpi = 300
)


# gene pathway network 


cnet <- cnetplot(
  ego_up,
  showCategory = 10
)

ggsave(
  "DR_Cimmune_gene_network.jpg",
  plot = cnet,
  width = 8,
  height = 8,
  dpi = 300
)

cnetplot(ego_up, showCategory = 8)


#====================================================================

# check sample distribution for each subgroup

table(meta.clean$subgroup)

# now do same things for control vs outbreak strain 


keep_con_outb <- meta.clean$subgroup %in% c('Outbreak TB strain','Control')

meta.clean_OB_C <- meta.clean[keep_con_outb,]

expr.OB_C <- expr[,keep_con_outb]
head(expr.OB_C)
head(meta.clean_OB_C)
# now convert the subgroup in meta.clean_OB_C into factor 

meta.clean_OB_C$subgroup <- factor(
  meta.clean_OB_C$subgroup,
  levels = c("Outbreak TB strain","Control")
)


# now run limma for the analysis of control vs Outbreak Strain 

library(limma)
# design matrix

design <- model.matrix(~subgroup, data = meta.clean_OB_C)

fit <- lmFit(expr.OB_C, design = design)
fit <- eBayes(fit)

# Get results
res_OB_C <- topTable(
  fit,
  coef = "subgroupOutbreak TB strain",
  number = Inf,
  adjust.method = "BH"
)

head(res)

table(head(res))
