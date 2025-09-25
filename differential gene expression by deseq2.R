# load the library
library('DESeq2')
library('parathyroidSE')

# data command load a dada object named parathyroidGenesSE 

data("parathyroidGenesSE")
# use assay() function to see the actual data 
# in this table each row represent ensemble gene id 
# each column represent sequenced rna library raw count

head(assay(parathyroidGenesSE))



# visualizing the metadata for the sample
colData(parathyroidGenesSE)


# number of technical replicates 

as.data.frame(colData(parathyroidGenesSE)[,c("sample","patient","treatment","time")])

# find the all the sample id
# calculate in which position same sample id exist
allColSamples<-colData(parathyroidGenesSE)$sample
sp<-split(seq(along=allColSamples),colData(parathyroidGenesSE)$sample)


countdata<-sapply(sp,function(columns)
  rowSums(assay(parathyroidGenesSE)[,columns,drop=FALSE]))
head(countdata)

# now manually perform the previous operation of the joining same
# count calue for the determined technical replicate 
a<-assay(parathyroidGenesSE)
countdata2<-cbind(a[,1:8],a[,9]+a[,10],a[,11],
                  a[,12]+a[,13],a[,14:22],a[,23]+a[,24],a[,25],a[,26]+a[,27])

# check coundata and countdata2 get the same result on not

all(countdata==countdata2)
print(countdata2)


# construct data object by rowdata, coldata, countdata

coldata<-colData(parathyroidGenesSE)[sapply(sp,`[`,1),]
rownames(coldata)<-coldata$sample
coldata

# keep the coldata column that actually need for analysis
coldata <- coldata[ ,c("patient","treatment","time")]
head(coldata)

# bring row metadata
rowdata <- rowData(parathyroidGenesSE)
rowdata

#now construct full dataset by DESeqDataSetFromMatrix()

ddsFull <- DESeqDataSetFromMatrix(
  countData = countdata2,
  colData = coldata,
  rowData = rowdata,
  design = ~ patient + treatment,
)

# now subset the data accordingly
# keep the data where treatment is control and dpn
# and time is 48 hour
dds <- ddsFull[ , colData(ddsFull)$treatment %in% c("Control","DPN") & 
                  colData(ddsFull)$time == "48h"]

# now selected data is refactored in dds
dds$patient <- factor(dds$patient)
dds$treatment <- factor(dds$treatment)


# relevel the treatment column with respect to Control
dds$treatment <- relevel(dds$treatment, "Control")


#we have to check the metadata of the dds dataset
colData(dds)

#run the dseq
dds <- DESeq(dds)



# result of the analysis
res <- results(dds)

#as it is a dataframe object, it has metadata
mcols(res)

# check for how many gene p value is less than 0.05
sum(res$pvalue<0.05, na.rm = TRUE)

sum(res$padj<0.1, na.rm = TRUE)

#subset the gene with adjp<0.1

resSort<-res[which(res$padj<0.1),]

# find out most downregulated gene from the adjusted pvalue sort
head(resSort[order(resSort$log2FoldChange),])

# visualization of the result with MAplot
plotMA(dds, ylim=c(-1.5,1.5))


plotDispEsts(dds)
sqrt(0.01)
