####  DESeq2 Script for Differential Gene Expression Analysis in 
# Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

## Use the Session menu to set working directory To Source File Directory

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
### How you inport this will depend on what your final output was from the mapper/counter that you used.
## this works with output from PrepDE.py from Ballgown folder.
countdata <- as.matrix(read.csv("/Users/waqartabish/Downloads/gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)
# Remove ".fastq" from column names if any
colnames(countdata) <- gsub(".fastq","",colnames(countdata))

# Clean the row names remove from the start till "|"
rownames(countdata) <- gsub(".*\\|","",rownames(countdata))
#### IF necessary,depending on what program made the count matrix. Remove the unwanted row (the * and zero row)
#countdata<- countdata[-1,c(-1)]
# OR   Remove the unwanted column (length)
#countdata<- countdata[,c(-1)]
#dim(countdata)
#head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <- read.table("/Users/waqartabish/Downloads/Proposed_Samples_Meta_Data.txt", 
                      header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE, row.names=1)
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
# which ones are not in the countdata?
rownames(coldata)[!rownames(coldata) %in% colnames(countdata)]

# only keep the samples in the coldata which are in the countdata
coldata <- coldata[colnames(countdata),]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~Size)
#look at it
dds



#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
## You can play around with this number to see how it affects your results!
dds <- dds[ rowSums(counts(dds)) > 10, ]
# look.  How many genes were filtered out?
dds

## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds$Size <- factor(dds$Size, levels=c("Small","Large"))


######     1.4 Differential expression analysis
### Question 2. Look at the manual - what is happening at this point?
dds <- DESeq(dds)
res <- results(dds)
res






###    1.5.1 MA-plot
##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
## Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down
plotMA(res, main="DESeq2", ylim=c(-8,8))

##  Write your results to a file 
write.csv(as.data.frame(resOrdered), file="Dog_DGESeq_results.csv")  

## 2.1.2 Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

### Heatmap of the count matrix
#library("genefilter")
# not taking into account the genes starting with "LOC"


topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

library("pheatmap")
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Size", "sex")])
df <- as.data.frame(colData(dds)[,c("Size","sex")])
pheatmap(mat, annotation_col = anno)


## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
#  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
#  nt <- normTransform(dds) # defaults to log2(x+1)
#  log2.norm.counts <- assay(nt)[select,]
#   df <- as.data.frame(colData(dds)[,c("treatment","type")])
#  pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#          cluster_cols=TRUE, annotation_col=df)

# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df) 


#2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Size)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# 2.2.3 Principal component plot of the samples
plotPCA(rld, intgroup=c("Size"))



############ Preparing Data for GSEA and Cytoscape.  #############

## Import the DGE results file make sure the gene model name is 'gene_id' to match annotation file
DGEresults <- read.csv("Dog_DGESeq_results.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)


############################# Make ranked list for GSEA ####################
## example aa <- within(resOrdered, z <- x + y - 2)
DGE_Anno_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Anno_Rank 
# first column name change to "Name" to match GSEA input
names(DGE_Anno_Rank)[1]<- "Name"
#subset the results so only Gene Name and rank
DGErank = subset(DGE_Anno_Rank, select = c(Name,rank) )
DGErank

#sebset the results so only Gene Name and rank
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)



### IF NEEDED....remove the "gene-" from row names
## https://stackoverflow.com/questions/39897155/replace-and-remove-part-of-string-in-rownames/39897315  "URS000075AF9C-snoRNA_GTATGTGTGGACAGCACTGAGACTGAGTCT"    to   "snoRNA"
## We can use gsub to match one of more characters that are not a - ([^-]+) from the start (^) of the string followed by 
## a - or (|) one or more characters that are not an underscore ([^_]+) until the end of the string ($)  and replace it with blanks ("").
## gsub("^[^-]+-|_[^_]+$", "", rownames(df))

#gene <-gsub("^[^-]+-", "", rownames(DGErank))
#DGErankIDs  <-cbind(gene,DGErank)
#head(DGErankIDs)
#summary(DGErankIDs)  

#write.csv(as.data.frame(DGErankIDs), file="DGErankIDs.csv", row.names=FALSE)  
write.table(as.data.frame(DGErank_withName), file="Dog_DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  

####  We also need the normalized expression DATA
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)
# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
gene <-gsub("^[^-]+-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(gene,NormTransExp)
head(NormTransExpIDs)

#write.csv(as.data.frame(NormTransExpIDs), file="NormTransExpressionData.csv", row.names=FALSE)  
write.table(as.data.frame(NormTransExpIDs), file="Dog_NormTransExpIDs.txt", quote=FALSE, row.names=FALSE, sep = "\t")  

