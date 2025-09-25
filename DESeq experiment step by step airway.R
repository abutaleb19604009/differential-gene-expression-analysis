# how to use differential gene expression analysis by deseq2

# load the dataset frist 

library(dplyr) 
library(tidyverse)
library(airway)

# get the data
data("airway")
airway
sample.info<- as.data.frame(colData(airway))

sample.info <- sample.info[,c(2,3)]

# rewrite the column content

sample.info$dex<- gsub('untrt','untreated',sample.info$dex)
sample.info$dex<-  gsub('trt','treated', sample.info$dex)

# rewrite the comumn name
names(sample.info) <- c('cellLine', 'dexamethasone')

# write metadata into a csv file 

write.table(sample.info, file = "sample.info.csv", sep = ',', col.names = T, row.names = T, quote = F)

# extract the count matrix from airway data object using assay function 
countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

# read the countsData, that has been written just 

countsData <- read.csv('counts_data.csv')

# read the metadata 
colData <- read.csv('sample.info.csv')


# check the col name in countmatrix present in rowname in coldata
all(colnames(countsData) %in% rownames(colData))

# check their order

all(colnames(countsData)==rownames(colData))



# construct the deseq data set as dds
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData= countsData,
  colData = colData,
  design = ~dexamethasone
)

# prefiltering in dds data before going for analysis by removing low gene expression 
# calculate row sum and keep only that row with row morethan 10

keep<- rowSums(counts(dds))>=10
dds <- dds[keep,]
counts(dds)

# set control and experiment level to deseq, factor level

dds$dexamethasone<- relevel(dds$dexamethasone, ref = 'untreated')
# if the dataset has techincal replicate, we need to collapse them before differential 
# gene expression analysis, so that each row contain only one gene, never collapse biological
# replicalte

# run deseq function 

dds <- DESeq(dds)

# find the results
res<- results(dds)
res

# summary of result
summary(res)

# see how many of gene up, down regulated based on p-value 
res0.01<- results(dds, alpha = 0.01)

summary(res0.01)


# visualizing the results with plot 
# MAplot
plotMA(dds)


