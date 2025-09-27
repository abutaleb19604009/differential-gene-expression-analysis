# differential gene expression analysis with edgeR


# set the working directory
setwd("~/Quantum biology/edgeR")

# download the gene wise count data from ncbi geo with GSE60450
# load the dataset in R studio 

GenewiseCounts<- read.delim('GSE60450_Lactation-GenewiseCounts.txt', row.names = 'EntrezGeneID')

colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)
head(GenewiseCounts)

targets <- read.delim("targets.txt.txt", stringsAsFactors=FALSE)
# construct DGEList with count matrix, sample info, and gene annotation info

#BiocManager::install('edgeR')
library(edgeR)
# Create the DGEList correctly
group <- paste(targets$CellType, targets$Status, sep=".")
group <- factor(group)
table(group)



y <- DGEList(counts = GenewiseCounts[,-1],
                           group = group,
                           genes = GenewiseCounts[,1, drop=FALSE])


options(digits = 3)

y$samples


# add gene annotation data
BiocManager::install('org.Mm.eg.db')
library(org.Mm.eg.db)

#crete a new column in y$gene list and add geneid  wiht column name SYMBOL
# value of the column as ENTREZID
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype = 'ENTREZID', column='SYMBOL')

# There are some na values in the count,  
#Entrez Ids that no longer have official gene symbols are dropped from the analysis
# keep all row where symbol value is not na and all column

y <- y[!is.na(y$genes$Symbol),]

dim(y)

#filter the genes that are expressed more than 0.5 count per million for at least 2 library

keep<- rowSums(cpm(y)>0.5)>=2
table(keep)


y <- y[keep,]



# calculate the normalization factor for each sample 

y <- calcNormFactors(y)

y$samples



plotMD(y, column=1)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)



y <- estimateDisp(y, design, robust=TRUE)

plotBCV(y)



