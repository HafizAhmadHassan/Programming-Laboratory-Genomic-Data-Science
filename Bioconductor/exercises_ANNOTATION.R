###### BSGENOME
library(BSgenome)
## EXERCISE 1
# - Show which genomes are available in the BSgenome package
# - Install and load the "BSgenome.Ecoli.NCBI.20080805" genome
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("BSgenome.Ecoli.NCBI.20080805")
# - Show informations about the genome
# - Show the names of the sequences of the genome
# - Put the first and the second sequence in a DNAStringSet "set"
# - Show a View with the first 100 elements of the first sequence of set
# - Show 100 random Views from the second sequence of set
#   - the width of each view must be 50
# - Show the base frequencies of the previous views 
# - Show the letter frequency of "ACG" in the sequences of set
# - Show the letter frequency of "ACG" in Ecoli
# - Search the sequence "ATGATAGATAG" in all sequences of Ecoli  


## EXERCISE 2
# Exploit the package BSgenome to retrieve the mus musculus genome (latest release) from UCSC.
# (a)For the chromosome with maximum length retrieve the CG content percentage.
# (b) Do the same for all chromosomes using bsapply.
# (c) Check if pattern d<-DNAString("AATC") occurs in a chromosome of your choice. 
# (d) Count the occurences of the reverse of the pattern in all chromosomes and save it in n_oc
# (e) Now consider the chromosome with the lowest number of occurrences (but greater than 0) and
# retrieve the sequence that starts at the start of the first occurence and which ends in the last position
# of the last occurrence.


###### ANNOTATIONHUB 

## EXERCISE 1
# Using the AnnotationHub package retrieve a gtf file for Homo Sapiens
#   - use the Ensembl database
#   - consider the GRCh37 genome (h19)
#   - use the makeTxDbFromGRanges function to convert the GRanges 
#     obtained from AnnotationHub to a TxDb object
#   - show the following features from the TxDb object:
#       - genes
#       - exons
#       - transcripts
#   - Select a gene from the TxDb object and find its transcripts
#   - show the exons grouped by genes


## EXERCISE 2
# Use the AnnotationHub to extract UCSC data that is from Homo sapiens and also specifically from the hg19 genome. What happens to the hub object as you filter data at each step?


## EXERCISE 3 
# Now that you have basically narrowed things down to the hg19 annotations from UCSC genome browser, lets get one of these annotations. Find the oreganno track and save it into a local variable.


## EXERCISE 4
# Load the org.Hs.eg.db package
# Look at the help page for the different columns and keytypes values with: help("SYMBOL"). 
# Now use this information to look up the entrez gene and chromosome for the gene symbol "MSX2" (from org.Hs.eg.db object).


## EXERCISE 5
# Extract all of the gene symbols and their associated entrez gene ids from the org.Hs.eg.db package. 
# Then check the symbols for redundancy.
# Hint: use duplicated function to subset the redundant symbols


## EXERCISE 6
# Use the accessors for the TxDb.Hsapiens.UCSC.hg19.knownGene package to retrieve the gene id,
# transcript name and transcript chromosome for all the transcripts. 
# Do this using both the select() method and also using the transcripts() method. 
# What is the difference in the output?


## EXERCISE 7
# Use the src_organism object with the transcripts method to look up the entrez gene IDs for all gene symbols that 
# contain the letter 'X'.


## EXERCISE 8
# Pull down GO terms for entrez gene id "1" from human by using the ensembl "hsapiens_gene_ensembl" dataset.

