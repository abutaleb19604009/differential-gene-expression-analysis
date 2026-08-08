# install packages and dependencies

install.packages('BiocManager')
BiocManager::install('Biobase')
library(Biobase)

#------------------------------------------------------------------------------
# open connection for data
con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = con)
# close the connection
close(con)

#-------------------------------------------------------------------------------

# take a short name
bm <- bodymap.eset

# extracting expression, phenodata, and feature data

exp_data <- exprs(bm)
peno_data <- pData(bm)
features_data <- fData(bm)

#-----------------------------------------------------------------------------
# check the dimention of the expression, phenodata, and features data
head(exp_data, n =5)
dim(exp_data)
# 52580    19

head(peno_data, n = 5)
dim(peno_data)
# 19  6

head(features_data, n =5)
dim(features_data)
# 52580     1

#------------------------------------------------------------------------------
# start with the exploratory data analysis

# improve visibility of scatter plot
install.packages("devtools")
devtools::install_github("alyssafrazee/RSkittleBrewer")
library(RSkittleBrewer)

trop <- RSkittleBrewer("tropical")
palette(trop)
par(pch = 19)

# load the required libraries 
BiocManager::install("org.Hs.eg.db")

library(gplots) 
library(devtools) 
library(Biobase)
library(RSkittleBrewer) 
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr) 

# start with tabulation of data, this helps to identify the distribution of 
# check the phenotypic imbalance among experimental group

table(peno_data$gender, peno_data$race)


# examine distribution of expression data for each sample

summary(exp_data)
# this revealed that genomic data are highly skewed 
#------------------------------------------------------------------------------

#check for the missing value in phenotypic data

table(peno_data$age, useNA='ifany')

table(is.na(peno_data$age))

# checking missing value in expression data

# along each row, gene in the expression data
missing_gene <- rowSums(is.na(exp_data))
table(missing_gene)

# along each column or along sample

missing_sample <- colSums(is.na(exp_data))

table(missing_sample)

#------------------------------------------------------------------------------

# now start the visualization 

# start with boxplot for expression data
# boxplot need a dataframe, first logtransform, then convert it to long format

expr_long <- as.data.frame(log2(exp_data)+1)
expr_long$Gene <- rownames(expr_long)

boxplot(log2(exp_data)+1, col=2, range=0)

# histogram for sample1
# with log transform
hist(log2(exp_data[ ,1]+1))
# without log transform
hist(exp_data[ ,1]+1)


# density plot for sample 1

plot(density(log2(exp_data[ ,1]+1)))
lines(density(log2(exp_data[ ,3]+1)))

# quantile plot

qqplot(log2(exp_data[ ,1]+1),
       log2(exp_data[ ,2]+1), col=3))
abline(c(0,1))

# MA Plot or Bland-Altman plot

mm <- (log2(exp_data[,1]+1)-log2(exp_data[,2]+1))
aa <- (log2(exp_data[,1]+1)+log2(exp_data[,2]+1))
plot(aa,mm, col =2)

#--------------------------------------------------------------------------
# filter low expressing gene

#only genes with an average expression greater than 1 are retained. Genes 
#with extremely low counts are removed from further visualization

exp_data <- as.data.frame(exp_data)
fil_data <- filter(exp_data, rowMeans(exp_data)>1)
#now visualize the exp data
boxplot(as.matrix(log2(fil_data + 1)),
        col = 2)

#------------------------------------------------------------------------------
# independent external validation 


# heatmap visualization 

