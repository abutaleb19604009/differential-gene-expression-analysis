# differential gene expression analysis practice from already 
# normalized data
# load the dataset frist 

library(dplyr) 
library(tidyverse)
library(GEOquery)

# read the normalized count data

data <- read.csv('GSE183947_fpkm.csv')
count.data <- data
# load the metadata from the GEO 
gse <- getGEO(GEO ='GSE183947', GSEMatrix = TRUE )
metadata <- pData(phenoData(gse[[1]]))

# select only required column from the metadata table 

metadata.subset <- select(metadata, c(1,10,11,17))


# rename of the column to easy analysis


metadata.subset <-metadata.subset %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastatis = characteristics_ch1.1)

# mutate the column value in easy format
metadata.subset<- metadata.subset %>%
  mutate(tissue= gsub("tissue: breast tumor","tumor",tissue)) %>%
  mutate(tissue= gsub("tissue: normal breast tissue","normal",tissue)) %>%
  mutate(metastatis= gsub("metastasis: yes","yes",metastatis)) %>%
  mutate(metastatis= gsub("metastasis: no","no",metastatis))
  


# reshape my count.data from wide format to long data
data.long<- count.data %>%
  rename(gene=X)%>%
  gather(key = 'Samples', value = 'FPKM',-gene)

# join long data and metadata.subset

joined.data<- data.long%>%
  left_join(.,metadata.subset, by=c('Samples'='description'))
  

# explore the joined data

joined.data%>%
  filter(gene=='BRCA1'| gene == 'BRCA2')%>%
  group_by(gene, tissue)%>%
  summarise(mean_fpkm = mean(FPKM),
            median_fpkm=median(FPKM))%>%
  head()



