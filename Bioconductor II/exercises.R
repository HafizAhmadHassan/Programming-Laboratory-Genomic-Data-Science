# Exercise 1

# Consider the series GSE87495 from the GEO platform
# This series contains 6 samples of human cells in two conditions:
#   - 3 samples treated with doxycycline
#   - 3 samples with no doxycycline
# Read the csv file "GSE87495_counts.csv" into a dataframe df
# The last 6 columns contain the counts for each sample  
colnames(df)
dim(df)
# Change the df column names related to the samples using the names in the 
# variable "samples" defined previously
df[,4:9]<-samples 
df
rownames(df)<- df[,3]
df2<-df[,3:9]

coldata<-data.frame(STATUS=samples)
rownames(df2)<-df2[,1]
m<-as.matrix(df2)

SummarizedExperiment(assays=list(counts=as.matrix(df2)), colData=coldata)

# Create SummarizedExperiment object SE: subset df (counts and genes) and assign colnames and rownames to the ones created.
# assay name must be "exp1"
exp1 <-df[3]
# colData is a 6 x 1 dataframe with a STATUS variable "no_dox" or "with_dox" based on colnames
coldata<-data.frame(STATUS=samples)
coldata
  # Add dataType="rseq" as metadata.
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(methods)
# Shows the structure and its elements
  
row.names(df)<-df[4]
names(df)
df <- data.frame(chr="chr2", start = 11:15, end = 12:16,
                 strand = c("+", "-", "+", "*", "."), expr0 = 3:7,
                 expr1 = 8:12, expr2 = 12:16,
                 row.names = paste0("GENE", letters[5:1]))
colnames(df)
row.names(df)<-df[,3]
exRSE <- makeSummarizedExperimentFromDataFrame(df)

exRSE

assay(exRSE)

rowRanges(exRSE)

nrows <- 200
ncols <- 6
counts1 <- matrix(runif(nrows * ncols, 1, 2000), nrows)
counts1
counts2 <- matrix(runif(nrows * ncols, 1, 50), nrows)

colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
colData

df = read.csv("GSE87495_counts.csv")
df
# We want to change the names of the samples of the dataset
# Define the variable "samples"
samples = c("no_dox_0", "with_dox_1", "no_dox_2",
            "with_dox_3", "no_dox_4", "with_dox_5")

colnames(df)[4:9]<-samples

df
df2<-data.frame(df[,4:9])
df2
colData<-DataFrame(STATUS = c(0,1,0,1,0,1),
                   row.names=samples)
df
rowData <-df[,1:3]

se=SummarizedExperiment(assays=list(exp1=df2),
                        rowData = rowData, colData=colData)
se

exp1<-assay(se)

assay(se)$dataType <-"rseq"
metadata(se)$dataType<-"rseq"
rowData(se)
#-------------------------------------------------------------------

#Exercise 2

# Create a new SummarizedExperiment DOX with only sample with dox from the 
# SummarizedExperiment of the previous Exercise

colData<-DataFrame(STATUS = c(1,1,1),
                   row.names=samples[c(2,4,6)])
se2=SummarizedExperiment(assays=list(exp1=df2[,c(2,4,6)]),
                        rowData = rowData, colData=colData)
se2

# Create a new random assay(matrix) and add it to DOX (values from normal distribution with mean = 0 and sd = 1)
# Hint: you can create a new SummarizedExperiment with 2 assays or replace the assays(DOX) list
# Hint: dimensions of the assays must be the same
nrow(df)
m<-matrix(rnorm(nrow(df)*3),nrow=nrow(df))
se2=SummarizedExperiment(assays=list(exp1=df2[,c(2,4,6)],exp2=m),
                         rowData = rowData, colData=colData)
se2

exp1[rowData(se2)$ENTREZID == "100288801",]
# Select data of the gene with ENTREZID = 100288801

#-------------------------------------------------------------------

# Exercise 3

# Build a RangedSummarizedExperiment from the SummarizedExperiment SE (Exercise 1)
# adding to SE the following GRanges:
# - seqnames
# - ranges
# - strand
# - type



#Show the class of SE


# Donwload from AnnotationHub gtf file for Homo Sapiens 
library(AnnotationHub)
hub <- AnnotationHub()
hub2<-hub[hub$dataprovider=="Ensembl"]
hub2 <- hub2[hub2$species=="Homo sapiens"]

# use the Ensembl database consider the GRCh37 genome (h19)
re<-hub2[['AH7558']]

# Filter the SE row according to the GRanges obtained 
#from GRanges obtained from AnnotationHub 

# How much rows have been removed?

#-------------------------------------------------------------------

# Exercise 4

# From the MultiAssayExperiment 

library(MultiAssayExperiment)
#miniACC select only assay RNASeq2GeneNorm and miRNASeqGene
data(miniACC)
se1<-miniACC[["RNASeq2GeneNorm"]]
se2 <-miniACC[["miRNASeqGene"]]
# This two assays have rows in common? Show the results
rownames(se1)
rownames(se2)
# This two assays have columns in common? Show the results
intersect(colnames(se1),colnames(se2))
# show the number of male and female in the MultiAssayExperiment
experiments(miniACC)[c(1,5)]


# Select only female sample 
#by creating a new MultiAssayExperiment called female
fse<-miniACC[, miniACC$gender=='female']

# From female add a new experiment as matrix with row and col equals to 
#RNAseq but random data
experiments(miniACC)[c(1,2)]

fse
miniACC2 <- c(miniACC, rnaseq = matrix(rnorm(nrow(assays(miniACC)$RNASeq2GeneNorm)*ncol(assays(miniACC)$RNASeq2GeneNorm)),nrow = nrow(assays(miniACC)$RNASeq2GeneNorm),ncol=ncol(assays(miniACC)$RNASeq2GeneNorm)), mapFrom=1)

#-------------------------------------------------------------------


#Exercise 5


# Try to extract the first SummarizedExperiment 
#from miniACC using [ e [[
e1<-experiments(miniACC)[[1]]
e2<-experiments(miniACC)[1]
# Are there differences in the obtained results?
assay(e)
miniACC[[1]]
# Is it possible to extract samples with vital_status = 0?
# Motivate the answers

#no sample have no meetadata
e2[[1]]
# Using the two methods 
#[ and [[ in combination it is possible to set correctly
miniACC[1]
# colData(...)=colData(...)
n#o
# Is it possible to extract samples with vital_status = 0? Try

assay(e2)
# Then try extracting the SummarizedExperiment using the getWithColData method (type ?getWithColData)
# Is it possible to extract samples with vital_status = 0?
# Motivate the answer

#no


#-------------------------------------------------------------------

#Exercise 6

# Download TGCA Data for disease UCEC (only assay RNASeqGene and miRNASeqGene)
library(curatedTCGAData)
library(TCGAutils)


curatedTCGAData(diseaseCode = "UCEC")

ucec2=curatedTCGAData(diseaseCode = "UCEC", assays = "miRNASeqGene*", dry.run = FALSE)

ucec1=curatedTCGAData(diseaseCode = "UCEC", assays = "RNASeqGene*", dry.run = FALSE)


# Extract the SummarizedExperiment RNASeqGene 
#from this MultiAssayExperiment 

# with all the methods seen in the previous exercise. 
#Make sure colData is not empty.
colData(ucec1)
# Hint: print ncol(assay) and nrow(colData).
# Hint: use the fuction anyReplicated (?anyReplicated) 
# and mergeReplicates (?mergeReplicates)
anyReplicated(ucec1)
mergeReplicates(ucec1)
# Try doing this exercise again without using the marging features (mergeReplicates)

# Hint: print colnames() and rownames().
colnames(ucec1)
rownames()
# Hint: converting the colnames names of the SummarizedExperiment assay into the barcode
# containing only the patient information it is possible to identify (and remove) 
# the duplicates with the classic summarized operations. 
# TCGAbarcode(..., participant = TRUE)
# https://bioc.ism.ac.jp/packages/release/bioc/vignettes/TCGAutils/inst/doc/TCGAutils.html

(xbarcode <- head(colnames(ucec1)[["UCEC_RNASeqGene-20160128"]], 4L))

TCGAbarcode(xbarcode, participant = TRUE, sample = TRUE, portion = FALSE)

barcodeToUUID(xbarcode)
#-------------------------------------------------------------------

#Exercise 7
# https://www.bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html

# Create a MultiAssayExperiment with:
empty <- MultiAssayExperiment()

# - a SummarizedExperiment with 2 assays containing counts of 
#6 samples and 50 genes, samples information include 
#   Condition (control or case, randomly assigned with probability 50%) 
#and samples are named sample1, sample2...
exprss1 <- matrix(rnorm(300), ncol = 6,nrow=50,
                  dimnames = list(sprintf("ENST00000%i", sample(288754:290000, 50)),
                                  c("sample1","sample2","sample3","sample4","sample5","sample6")))
exprss2 <- matrix(rnorm(300), ncol = 6,
                  dimnames = list(sprintf("ENST00000%i", sample(288754:290000, 50)),
                                  c("sample1","sample2","sample3","sample4","sample5","sample6")))

exprss2


doubleExp <- list("Expr1"  = exprss1, "expr2" = exprss2)
simpleMultiAssay <- MultiAssayExperiment(experiments=doubleExp)
simpleMultiAssay

coll.data <- data.frame(Condition=sample(c("control", "case"),size=6,replace = TRUE),
                           
                           row.names=c("sample1","sample2","sample3","sample4","sample5","sample6"))

simpleMultiAssay2 <- MultiAssayExperiment(experiments=doubleExp,
                                          colData=coll.data)

# - a matrix with 30 rows (named with "hsa_" and a random number) 
#and 4 columns (called array1, array2...)
m <- matrix(rnorm(120), ncol = 4,
                  dimnames = list(sprintf("hsa_%i", sample(288754:290000, 30)),
                                  c("array1", "array2", "array3", "array4")))

m

# Create also a data.frame with patients data (sex and name) and 
patient.data <- data.frame(sex=c("M", "F", "M", "F"),
                           age=38:41,
                           row.names=c("array1", "array2", "array3", "array4"))
patient.data

#define a mapping for the SummarizedExperiment and matrix just created.

exprdat<-SummarizedExperiment(assays=list(counts=as.matrix(m)), colData=patient.data)
(exprmap <- data.frame(primary=rownames(patient.data)[c(1, 2, 4, 3)],
                      colname=c("array1", "array2", "array3", "array4"),
                      stringsAsFactors = FALSE))
listmap<-list(exprmap)
names(listmap)<-c("Mymapping")
dfmap<- listToMap(listmap)

myMultiAssay <- MultiAssayExperiment(list(ex1=exprss1,ex2=exprss2), coll.data, dfmap)

