#Exercise 1

library(EnrichmentBrowser)

# Consider the series GSE87495 from the GEO platform
# This series contains 6 samples of human cells in two conditions:
#   - 3 samples treated with doxycycline
#   - 3 samples with no doxycycline
# Read the csv file "GSE87495_counts.csv" into a dataframe df
# The last 6 columns contain the counts for each sample  

df <- read.csv("GSE87495_counts.csv")

# We want to change the names of the samples of the dataset
# Define the variable "samples"
samples = c("no_dox_1", "with_dox_2", "no_dox_3",
            "with_dox_4", "no_dox_5", "with_dox_6")
samples

# Change the df column names related to the samples using the names in the 
# variable "samples" defined previously

colnames(df)[4:9] <- samples
colnames(df)[4:9]

#Create SummarizedExperiment object: subset df (counts and genes) and assign colnames and rownames to the ones created.
# Hint: counts must be of class matrix, while the other can be dataframe.
# Hint: create a dataframe colData with a GROUP variable: set no_dox GROUP=0, with_dox=1

genes <- df[,1:3]
genes
counts <- as.matrix(df[,4:9])
counts

rownames(counts) = df[,1]

colnames(counts) = colnames(df[4:9])

rownames(genes)= df[,1]


#colData: rownames as colnames(counts) ;set no_dox GROUP=0, with_dox=1

colData = DataFrame(name=colnames(counts),status=substr(colnames(counts),1,nchar(colnames(counts))-2) ,row.names=colnames(counts))
colData$GROUP= ifelse(colData$status=="no_dox",0,1)

se = SummarizedExperiment(assays=counts, rowData=genes, colData=colData, metadata=list(dataType="rseq", annotation = "hsa"))
se

#normalize data ("normalize" function)

dd <- normalize(se)
dd

#Compute differentially expressed genes analysis (deAna function).
#hint: show rowData(se) before and after the calculation.

dd <-deAna(dd)

#set-bases enrichment
#1. create gene sets:
# 1.1 use the function getGenesets to download all KEGG pathways for Homo sapiens as gene sets.
# 1.2 use the function getGenesets to retrieve GO terms of a biological process (BP) ontology as defined in the GO.db annotation package.
#2. compute enrichment for both gene sets using a method of your choice
#3. show results (gsRanking and eaBrowse)

library("ReportingTools")

kegg.gs = getGenesets(org="hsa", db="kegg")
go.gs = getGenesets(org="hsa", db="go", go.onto="BP", go.mode="GO.db")

kegg.se = sbea(method = "ora", se=dd, gs=kegg.gs, perm = 0)
go.se = sbea(method = "ora", se=dd, gs=go.gs, perm = 0)

gsRanking(kegg.se)
gsRanking(go.se)
eaBrowse(kegg.se,nr.show=5)
eaBrowse(go.se,nr.show=5) 

#network-based enrichment
#1. compile a gene regulatory network from regulations in KEGG pathway databases (organism Homo Sapiens)
#2. compute network-based enrichment using a method of your choice
#3. show results (gsRanking) 
#4. plot Gene graph enrichment analysis (GGEA) for pathway "hsa04913_Ovarian_steroidogenesis" (with legend)

hsa.grn = compileGRN(org="hsa", db="kegg")
kegg.seGRN = nbea(method="ggea", se=dd, gs=kegg.gs, grn=hsa.grn)
gsRanking(kegg.seGRN)

ggeaGraph(gs=kegg.gs[["hsa04913_Ovarian_steroidogenesis"]], grn=hsa.grn, se=dd)
ggeaGraphLegend()

# combining results. calculate and display (eaBrowse) the results obtained by combining:
#1. set-based and network-based calculated first using KEGG gene sets
#2. two different set-based enrichment methods (gene sets from GO)
#3. two different methods of network-based enrichment

res.list = list(kegg.se, kegg.seGRN)
comb.res = combResults(res.list)
eaBrowse(comb.res, graph.view=hsa.grn, nr.show=5)

go.se2 = sbea(method = "safe", se=dd, gs=go.gs, perm = 0)
res.list2 = list(go.se, go.se2)
comb.res2 = combResults(res.list2)
eaBrowse(comb.res2, graph.view=hsa.grn, nr.show=5)

#---------------------------------------------------------------------------------------------------------------------

# Exerecise 2

# What are the pathways affected by the development of COVID-19?
# The file SARS-CoV-2_infection.csv contains data about the transcriptional response
# of human lung epithelial cells to SARS-CoV-2 infection.
# It contains 6 samples: 3 of them mock treated (control) and 3 of them infected with SARS-Cov-2
# Firstly, read the file into a dataframe df

df <- read.delim("SARS-CoV-2_infection.csv")


# Create a SummarizedExperiment from df

samples <- c("control_1", "control_2", "control_3",
             "infected_4", "infected_5", "infected_6")

 colnames(df)[2:7] <- samples

genes <- df[,1]
counts <- as.matrix(df[,2:7]) 
rownames(counts) = df[,1]
colnames(counts) = colnames(df[2:7])

colData = DataFrame(name=colnames(counts),status=substr(colnames(counts),1,nchar(colnames(counts))-2) ,row.names=colnames(counts))
colData$GROUP = c(rep(0,3), rep(1,3))

se <- SummarizedExperiment(assays = counts, rowData = genes, colData = colData)
se

# What's the current gene identifier?
# In order to perform enrichment later, we need the gene identifiers to be ENTREZID. Map them.

SE = idMap(se, org="hsa", from="SYMBOL", to="ENTREZID")

# Let's define the col data by creating a dataframe that has:
# - the name of the samples as row names
# - a column GROUP with three 0s and three 1s


# Normalize the data and perform differential expression analysis

dd <- normalize(SE)
dd <- deAna(dd)

# Use KEGG to perform a set-based enrichment analysis using ora and gsea

library("ReportingTools")

kegg.gs = getGenesets(org="hsa", db="kegg")
kegg.se = sbea(method = "ora", se=dd, gs=kegg.gs, perm = 0)
kegg.se2 = sbea(method = "gsea", se=dd, gs=kegg.gs, perm = 50)
gsRanking(kegg.se)
gsRanking(kegg.se2) 

# Use KEGG to perform a network-based enrichment analysis using spia

hsa.grn = compileGRN(org="hsa", db="kegg")
kegg.seGRN = nbea(method="spia", se=dd, gs=kegg.gs, grn=hsa.grn)
gsRanking(kegg.seGRN)

# Show results from both methods
# Combine the results obtained from the two methods and show the first 5 pathways
res.list = list(kegg.se, kegg.se2, kegg.seGRN)
comb.res = combResults(res.list)
eaBrowse(comb.res, graph.view=hsa.grn, nr.show=5)
#---------------------------------------------------------------------------------------------------------------------

#Exercise 3

# Download TGCA Data for disease UCEC (only assay RNASeqGene)

library(SummarizedExperiment)
library(MultiAssayExperiment)
library(curatedTCGAData)
library(TCGAutils)


UCEC.data <- curatedTCGAData(diseaseCode = "UCEC", assays = "RNASeqGene", dry.run = FALSE)

NameClin <- getClinicalNames("UCEC")

repli <- replicated(UCEC.data)
UCEC.data <- mergeReplicates(UCEC.data, replicates = repli)
colData(UCEC.data) <- colData(UCEC.data)[, 1:2]
vitalStat <- c("GROUP") #Corretto?
colnames(colData(UCEC.data))[2] <- vitalStat

#colData$GROUP <- colData$vital_status.x
se <- getWithColData(UCEC.data, 1L)
 
# Extract a SummarizedExperiment from MultiAssayExperiment just dowloaded using as colData only relevant columns 
# (use getClinicalNames, GROUP as vital status) and considering only RNASeqGene experiments 
# Hint: check for duplicates (replicated mergeReplicate); getWitColData to extract SE.

#Compute differentially expressed genes analysis (deAna function)

SE = idMap(se, org="hsa", from="SYMBOL", to="ENTREZID")
dd <- normalize(SE)
dd <- deAna(dd)

#compute set-based (only KEGG) and network-based enrichment and show results
# Hint: genes are stored as SYMBOL, it is necessary to convert them to ENTREZID before doing the analysis

kegg.gs = getGenesets(org="hsa", db="kegg")
kegg.se = sbea(method = "ora", se=dd, gs=kegg.gs, perm = 0)
gsRanking(kegg.se)
eaBrowse(kegg.se,nr.show=5)

hsa.grn = compileGRN(org="hsa", db="kegg")
kegg.seGRN = nbea(method="ggea", se=dd, gs=kegg.gs, grn=hsa.grn)
gsRanking(kegg.seGRN) 

