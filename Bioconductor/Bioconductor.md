---
title: " Bioconductor"
author: "Rosalba Giugno"
date: "11/11/2019"
output:
  html_document:
    keep_md: yes
    number_sections: no
    toc: yes
  slidy_presentation: default
---




# What is Bioconductor

Bioconductor provides tools for the analysis and comprehension of high-throughput genomic data. Bioconductor uses the R statistical programming language, and is open source and open development. It has two releases each year, and an active user community. Bioconductor is also available as an AMI (Amazon Machine Image) and a series of Docker images.

![Bioconductor](Bioconductor.png)

# Biology Pills

DNA, Chromosomes, and Genes  (4':39'')
https://www.youtube.com/watch?v=-i1_JagCL1U

Genome Sequencing 
https://www.youtube.com/watch?v=2JUu1WqidC4


Gene regulations
https://www.khanacademy.org/test-prep/mcat/biomolecules/gene-control/v/regulation-of-transcription

Post-trascriptional regulation
https://www.khanacademy.org/test-prep/mcat/biomolecules/gene-control/v/post-translational-regulation


# Install BiocManager


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

# Resources
https://www.bioconductor.org/install/

Packages

Help on line

Workflow

Documentation

# Install Packages 


```r
BiocManager::install(c("NAME PACKAGE"))
```

# Basic BioConductor Data Structures --IntervalRanges - IRanges

IRanges is a vector, that contains integer intervals. Sounds a little funny, but let's take an example. We construct an IRanges by using the IRanges constructor function, and we give two out of three arguments, start, end, and with. We only need two, because if we know two of them, the last one can be inferred. 


```r
BiocManager::install(c("IRanges"))
```




```r
library(IRanges)
```


```r
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         3         3
##   [2]         3         5         3
##   [3]         5         7         3
```
So here we have a start, an end, a width, and we can see the width column has been filled by knowing the start and the end. Here we construct another IRange by specifying the start and the width, and we get exactly the same object out. 


```r
ir2 <- IRanges(start = c(1,3,5), width = 3)
ir2
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         3         3
##   [2]         3         5         3
##   [3]         5         7         3
```

```r
all.equal(ir1, ir2)
```

```
## [1] TRUE
```

An IRanges consist of separate intervals; each interval is called a range. So ir1 above contains 3 ranges.

Methods: start(), end(), width() and also replacement methods.


```r
start(ir1)
```

```
## [1] 1 3 5
```
## Resize 

```r
width(ir2) <- 1
ir2
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         1         1
##   [2]         3         3         1
##   [3]         5         5         1
```
## Name

```r
names(ir1) <- paste("A", 1:3, sep = "")
ir1
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##          start       end     width
##      <integer> <integer> <integer>
##   A1         1         3         3
##   A2         3         5         3
##   A3         5         7         3
```
IRanges  have a single dimension, there are vectors no matrices

```r
dim(ir1)
```

```
## NULL
```

```r
length(ir1)
```

```
## [1] 3
```

## Subsetting 
Subsetting works like a vector

```r
ir1[1]
```

```
## IRanges object with 1 range and 0 metadata columns:
##          start       end     width
##      <integer> <integer> <integer>
##   A1         1         3         3
```

```r
ir1["A1"]
```

```
## IRanges object with 1 range and 0 metadata columns:
##          start       end     width
##      <integer> <integer> <integer>
##   A1         1         3         3
```
## Concatenate

Like vectors, you can concatenate two IRanges with the c() function

```r
c(ir1, ir2)
```

```
## IRanges object with 6 ranges and 0 metadata columns:
##          start       end     width
##      <integer> <integer> <integer>
##   A1         1         3         3
##   A2         3         5         3
##   A3         5         7         3
##              1         1         1
##              3         3         1
##              5         5         1
```



## Normal IRanges
A normal IRanges is a minimal representation of the IRanges viewed as a set. Each integer only occur in a single range and there are as few ranges as possible. In addition, it is ordered. Many functions produce a normal IRanges. Created by reduce().

Given the following IRange

```r
ir <- IRanges(start= c(1,3,7,9), end=c(4,4,8,10))
plotRanges(ir)
```

![](Bioconductor_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


```r
ir1 <- reduce(ir)
plotRanges(ir1) 
```

![](Bioconductor_files/figure-html/unnamed-chunk-16-1.png)<!-- -->
  
From some perspective, disjoin() is the opposite of reduce(). 

```r
ir2<- disjoin(ir)
plotRanges(ir2)
```

![](Bioconductor_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

## Manipulating IRanges, intra-ranges

Intra-range manipulations are manipulations where each original range gets mapped to a new range. Examples of these are: shift(), narrow(), flank(), resize(), restrict().

For example, resize() can be extremely useful. It has a fix argument controlling where the resizing occurs from.

Use fix="center" to resize around the center of the ranges.

```r
ira<- resize(ir, width = 1, fix = "start")
plotRanges(ira)
```

![](Bioconductor_files/figure-html/unnamed-chunk-18-1.png)<!-- -->


```r
ira<- resize(ir, width = 1, fix = "center")
plotRanges(ira)
```

![](Bioconductor_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

## Manipulating IRanges, as sets

Manipulating IRanges as sets means that we view each IRanges as a set  of integers; individual integers is either contained in one or more ranges or they are not. 
This is equivalent to calling reduce() on the IRanges first.

Once this is done, we can use standard:  union(), intersect(), setdiff(), gaps() 
between two IRanges (which all returns normalized  IRanges).


```r
par(mfrow = c(3,1))
ir1 <- IRanges(start = c(1, 3, 5), width = 1)
plotRanges(ir1)
ir2 <- IRanges(start = c(4, 5, 6), width = 1)
plotRanges(ir2)
ir3<- union(ir1, ir2)
plotRanges(ir3)
```

![](Bioconductor_files/figure-html/unnamed-chunk-20-1.png)<!-- -->


```r
ir4 <-intersect(ir1, ir2)
plotRanges(ir4)
```

![](Bioconductor_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

Because they return normalized IRanges, an alternative to union() is

```r
ir5 <- reduce(c(ir1, ir2))
plotRanges(ir5)
```

![](Bioconductor_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

## Overlaps
Finding (pairwise) overlaps between two IRanges is done by findOverlaps(). This function is very important and amazingly fast!

```r
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir1
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         3         3
##   [2]         4         7         4
##   [3]         8        10         3
```

```r
ir2 <- IRanges(start = c(3,4), width = 3)
ir2
```

```
## IRanges object with 2 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         3         5         3
##   [2]         4         6         3
```

```r
ov <- findOverlaps(ir1, ir2)
ov
```

```
## Hits object with 3 hits and 0 metadata columns:
##       queryHits subjectHits
##       <integer>   <integer>
##   [1]         1           1
##   [2]         2           1
##   [3]         2           2
##   -------
##   queryLength: 3 / subjectLength: 2
```

It returns a Hits object which describes the relationship between the two IRanges.
This object is basically a two-column matrix of indicies into the two IRanges.
The two columns of the hits object can be accessed by queryHits() and subjectHits()
(often used with unique()).

For example, the first row of the matrix describes that the first range of ir1 
overlaps with the first range of ir2. Or said differently, 
they have a non-empty intersection:

```r
intersect(ir2[subjectHits(ov)[1]],
            ir1[queryHits(ov)[1]])
```

```
## IRanges object with 1 range and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         3         3         1
```


The elements of unique(queryHits)  gives you the indices of the query ranges which actually had an overlap; you need unique because a query range may overlap multiple subject ranges.

```r
queryHits(ov)
```

```
## [1] 1 2 2
```

```r
unique(queryHits(ov))
```

```
## [1] 1 2
```

The list of arguments to findOverlaps() is long; there are a few hidden treasures here. For example, you can ask to only get an overlap if two ranges overlap by a certain number of bases.


```r
args(findOverlaps)
```

```
## function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any", 
##     "start", "end", "within", "equal"), select = c("all", "first", 
##     "last", "arbitrary"), ...) 
## NULL
```

## CountOverlaps
For efficiency, there is also countOverlaps(), which just returns the number of overlaps.
This function is faster and takes up less memory because it does not have to keep 
track of which ranges overlap, just the number of overlaps.

```r
countOverlaps(ir1, ir2)
```

```
## [1] 1 2 0
```
## Finding nearest IRanges

Sometimes you have two sets of IRanges and you need to know which ones are closest to each other. Functions for this include  nearest(), precede(), follow(). 


```r
ir1
```

```
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         3         3
##   [2]         4         7         4
##   [3]         8        10         3
```

```r
ir2
```

```
## IRanges object with 2 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         3         5         3
##   [2]         4         6         3
```

```r
nearest(ir1, ir2)
```

```
## [1] 1 1 2
```

# Basic BioConductor Data Structures-- GenomicRanges - GRanges




```r
BiocManager::install(c("GenomicRanges"))
library(GenomicRanges)
```

GRanges are like IRanges with strand and chromosome. 
Strand can be +, - and *. 
The value * indicates ‘unknown strand’ or ‘unstranded’. 
This value usually gets treated as a third strand, which is sometimes confusing to users (examples below).

A GRange is very similar to an IRanges with some additional having to do with chromosome 
and strength. Cromosones in GRanges are called seqnames. They get created with the GRanges constructor:


```r
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
                ranges = IRanges(start = c(1,3,5), width = 3))
gr
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr1       3-5      -
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

Natural accessor functions: strand(), seqnames(), ranges(), start(), end(), width().

Because they have strand, we now have operations which are relative 
to the direction of transcription (upstream(), downstream()):

## Intra-ranges-operations

Flank generates flanking ranges for each range in x. 
If start is TRUE for a given range, the flanking occurs at the start, otherwise the end.
The widths of the flanks are given by the width parameter.
The widths can be negative, in which case the flanking region is reversed so that 
it represents a prefix or suffix of the range in x. The flank operation is illustrated below 
for a call of the form flank(x, 3, TRUE), where x indicates a range in x and - 
indicates the resulting flanking region

```r
  flank(gr, 2, start = TRUE)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1      -1-0      +
##   [2]     chr1       6-7      -
##   [3]     chr1       3-4      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
Other operations are 

```r
shift(gr, 5)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       6-8      +
##   [2]     chr1      8-10      -
##   [3]     chr1     10-12      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
resize(gr, 30)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1      1-30      +
##   [2]     chr1     -24-5      -
##   [3]     chr1      5-34      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

to get help ?"intra-range-methods" 

## Seqinfo
GRanges operates within a universe of sequences (chromosomes/contigs) and their lengths.
This is described through seqinfo:

```r
seqinfo(gr)
```

```
## Seqinfo object with 1 sequence from an unspecified genome; no seqlengths:
##   seqnames seqlengths isCircular genome
##   chr1             NA         NA   <NA>
```

```r
seqlengths(gr) <- c("chr1" = 10)
seqinfo(gr)
```

```
## Seqinfo object with 1 sequence from an unspecified genome:
##   seqnames seqlengths isCircular genome
##   chr1             10         NA   <NA>
```

```r
seqlevels(gr)
```

```
## [1] "chr1"
```

```r
seqlengths(gr)
```

```
## chr1 
##   10
```
Especially the length of the chromosomes are used by some functions. 
For example gaps() return the stretches of the genome not covered by the GRanges.

```r
gaps(gr)
```

```
## GRanges object with 5 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1         4      +
##   [2]     chr1      8-10      +
##   [3]     chr1       1-2      -
##   [4]     chr1      6-10      -
##   [5]     chr1      1-10      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome
```
In this example, we know that the last gap stops at 10, because that is the length of the chromosome. Note how a range on the * strand appears in the result.

Let us expand the GRanges with another chromosome

```r
seqlevels(gr) <- c("chr1", "chr2")
seqnames(gr) <- c("chr1", "chr2", "chr1")
gr
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr2       3-5      -
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 2 sequences from an unspecified genome
```
When you sort() a GRanges, the sorting order of the chromosomes is determined by 
their order in seqlevel. This is nice if you want the sorting “chr1”, “chr2”, …, “chr10”, …

```r
sort(gr)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr1       5-7      +
##   [3]     chr2       3-5      -
##   -------
##   seqinfo: 2 sequences from an unspecified genome
```

The sort is sensible to the order of defined seqlevels

```r
seqlevels(gr) <- c("chr2", "chr1")
sort(gr)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr2       3-5      -
##   [2]     chr1       1-3      +
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 2 sequences from an unspecified genome
```


You can associate a genome with a GRanges.

```r
genome(gr) <- "hg19"
gr
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr2       3-5      -
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 2 sequences from hg19 genome
```

```r
seqinfo(gr)
```

```
## Seqinfo object with 2 sequences from hg19 genome:
##   seqnames seqlengths isCircular genome
##   chr2             NA         NA   hg19
##   chr1             10         NA   hg19
```

This becomes valuable when you deal with data from different genome versions (as we all do),
because it allows R to throw an error when you compare two GRanges from different genomes, like

```r
gr2 <- gr
genome(gr2) <- "hg18"
findOverlaps(gr, gr2)

ERROR in mergeNameAtomicVector
```

## Basic GRanges Usage -- DataFrame
DataFrame class is very similar to the base data.frame class from R, but it allows columns of any class, provided a number of required methods are supported. For example, DataFrame can have IRanges as columns, unlike data.frame:


```r
ir <- IRanges(start = 1:2, width = 3)
df1 <- DataFrame(iranges = ir)
df1
```

```
## DataFrame with 2 rows and 1 column
##     iranges
##   <IRanges>
## 1       1-3
## 2       2-4
```

```r
df1$iranges
```

```
## IRanges object with 2 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         3         3
##   [2]         2         4         3
```
In the data.frame case, the IRanges gives rise to 4 columns, whereas it is a single column when a DataFrame is used.

```r
df2 <- data.frame(iranges = ir)
df2
```

```
##   iranges.start iranges.end iranges.width
## 1             1           3             3
## 2             2           4             3
```


##  Metadata

GRanges (unlike IRanges) may have associated metadata. This is immensely useful. 
The formal way to access and set this metadata is through values or elementMetadata or mcols, like

```r
gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr1       3-5      -
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))
gr
```

```
## GRanges object with 3 ranges and 1 metadata column:
##       seqnames    ranges strand |     score
##          <Rle> <IRanges>  <Rle> | <numeric>
##   [1]     chr1       1-3      + |       0.1
##   [2]     chr1       3-5      - |       0.5
##   [3]     chr1       5-7      + |       0.3
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

A much easier way to set and access metadata is through the $ operator to access

```r
gr$score
```

```
## [1] 0.1 0.5 0.3
```
to add

```r
gr$score2 = gr$score * 0.2
gr
```

```
## GRanges object with 3 ranges and 2 metadata columns:
##       seqnames    ranges strand |     score    score2
##          <Rle> <IRanges>  <Rle> | <numeric> <numeric>
##   [1]     chr1       1-3      + |       0.1      0.02
##   [2]     chr1       3-5      - |       0.5       0.1
##   [3]     chr1       5-7      + |       0.3      0.06
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## FindOverlaps

findOverlaps works exactly as for IRanges. 

```r
gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*",
               ranges = IRanges(start = c(1, 3, 5), width = 3))
gr2
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      *
##   [2]     chr2       3-5      *
##   [3]     chr1       5-7      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

```r
gr
```

```
## GRanges object with 3 ranges and 2 metadata columns:
##       seqnames    ranges strand |     score    score2
##          <Rle> <IRanges>  <Rle> | <numeric> <numeric>
##   [1]     chr1       1-3      + |       0.1      0.02
##   [2]     chr1       3-5      - |       0.5       0.1
##   [3]     chr1       5-7      + |       0.3      0.06
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
findOverlaps(gr, gr2)
```

```
## Hits object with 4 hits and 0 metadata columns:
##       queryHits subjectHits
##       <integer>   <integer>
##   [1]         1           1
##   [2]         2           1
##   [3]         2           3
##   [4]         3           3
##   -------
##   queryLength: 3 / subjectLength: 3
```


## SubsetByOverlaps
A common operation is to select only certain ranges from a GRanges which overlap something else.  Enter the convenience function  subsetByOverlaps

```r
subsetByOverlaps(gr, gr2)
```

```
## GRanges object with 3 ranges and 2 metadata columns:
##       seqnames    ranges strand |     score    score2
##          <Rle> <IRanges>  <Rle> | <numeric> <numeric>
##   [1]     chr1       1-3      + |       0.1      0.02
##   [2]     chr1       3-5      - |       0.5       0.1
##   [3]     chr1       5-7      + |       0.3      0.06
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
## MakeGRanges Fromn DataFrame

A common situation is that you have data which looks like a GRanges but is really stored 
as a classic data.frame, with chr, start etc. The makeGRangesFromDataFrame converts this data.frame into a GRanges. An argument tells you whether you want to keep any additional columns.

```r
df <- data.frame(chr = "chr1", start = 1:3, end = 4:6, score = 7:9)
df
```

```
##    chr start end score
## 1 chr1     1   4     7
## 2 chr1     2   5     8
## 3 chr1     3   6     9
```

```r
makeGRangesFromDataFrame(df)
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-4      *
##   [2]     chr1       2-5      *
##   [3]     chr1       3-6      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
```

```
## GRanges object with 3 ranges and 1 metadata column:
##       seqnames    ranges strand |     score
##          <Rle> <IRanges>  <Rle> | <integer>
##   [1]     chr1       1-4      * |         7
##   [2]     chr1       2-5      * |         8
##   [3]     chr1       3-6      * |         9
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Seqinfo 

```r
BiocManager::install(c("GenomeInfoDb"))
```

```
## Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.3 (2020-02-29)
```

```
## Installing package(s) 'GenomeInfoDb'
```

```
## Warning: package 'GenomeInfoDb' is in use and will not be installed
```

```
## Installation path not writeable, unable to update packages: boot, class,
##   KernSmooth, MASS, nlme, nnet, spatial, survival
```

```
## Old packages: 'purrr', 'RCurl', 'tibble'
```

```r
library(GenomeInfoDb)
```
## Drop and keep seqlevels
It is common to want to remove seqlevels from a GRanges object. 
Use  pruning.mode="coarse", the default value is  pruning.mode="error" meaning it returns an error and does not drop levels

```r
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
gr
```

```
## GRanges object with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-4      *
##   [2]     chr2       2-5      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

```r
seqlevels(gr, pruning.mode="coarse") <- c("chr1")
gr
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-4      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

You can also just get rid of weird looking chromosome names with keepStandardChromosomes().

```r
gr <- GRanges(seqnames = c("chr1", "c5hrU34"),
              ranges = IRanges(start = 1:2, end = 4:5))
gr
```

```
## GRanges object with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-4      *
##   [2]  c5hrU34       2-5      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

```r
keepStandardChromosomes(gr,pruning.mode="coarse")
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-4      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Changing style
It is an inconvenient truth that different online resources uses different naming convention 
for chromosomes. This can even be different from organism to organism. 
For example, for the fruitfly (Drosophila Melanogaster) NCBI and Ensembl
uses “2L” and UCSC uses “chr2L”. But NCBI and Ensembl differs on some contigs: NCBI uses “Un” and Ensembl used “U”.

```r
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 2))
gr
```

```
## GRanges object with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-2      *
##   [2]     chr1       2-3      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
#Let us remap
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
gr <- renameSeqlevels(gr, newStyle)
gr
```

```
## GRanges object with 2 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]        1       1-2      *
##   [2]        1       2-3      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

This can in principle go wrong, if the original set of seqlevels are inconsistent
(not a single style).
The GenomeInfoDb also contains information for dropping / keeping various classes of chromosomes:
  

## Rle (run-length-encoded) vector

We will discuss a data representation class called Rle (run length encoding). 
This class is great for representation genome-wide sequence coverage.
In high-throughput sequencing, coverage is the number of reads overlapping each base.  In other words, it associates a number (the number of reads) to every base in the genome.

This is a fundamental quantity for many high-throughout sequencing analyses.  For variant calling (DNA sequencing) it tells you how much power (information)  you have to call a variant at a given location.  For ChIP sequencing it is the primary signal;  areas with high coverage are thought to be enriched for a given protein. A file format which is often used to represent coverage data is Wig or the modern version BigWig.


Look at the name of the columns og GRanges. Look at Rle.  

```r
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
gr
```

```
## GRanges object with 10 ranges and 0 metadata columns:
##        seqnames    ranges strand
##           <Rle> <IRanges>  <Rle>
##    [1]     chr1       1-3      *
##    [2]     chr1       2-4      *
##    [3]     chr1       3-5      *
##    [4]     chr1       4-6      *
##    [5]     chr1       5-7      *
##    [6]     chr1       6-8      *
##    [7]     chr1       7-9      *
##    [8]     chr1      8-10      *
##    [9]     chr1      9-11      *
##   [10]     chr1     10-12      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```


An Rle (run-length-encoded) vector is a specific representation of a vector. The IRanges package implements support for this class. Watch out: there is also a base R class called rle which has much less functionality.

The run-length-encoded representation of a vector, represents the vector as a set of distinct runs with their own value. Let us take an example


```r
rl <- Rle(c(1,1,1,1,2,2,3,3,2,2))
rl
```

```
## numeric-Rle of length 10 with 4 runs
##   Lengths: 4 2 2 2
##   Values : 1 2 3 2
```

```r
runLength(rl)
```

```
## [1] 4 2 2 2
```

```r
runValue(rl)
```

```
## [1] 1 2 3 2
```

```r
as.numeric(rl)
```

```
##  [1] 1 1 1 1 2 2 3 3 2 2
```

This is a very efficient representation if the vector is very long there are a lot of consecutive elements with the same value. This is especially useful for genomic data which is either piecewise constant, or where most of the genome is not covered (eg. RNA sequencing in mammals).

In many ways Rles function as normal vectors, you can do arithmetic with them, transform them etc. using standard R functions like + and  log2.

There are also RleList which is a list of Rles.  This class is used to represent a genome wide coverage track where each element of the list is a different chromosome.


### Useful functions for Rle

A standard usecase is that you have a number of regions (say IRanges) and you want to do something to your Rle over each of these regions. Enter aggregate().
Compute the mean in the regions of the rle vector (that is a compressed vector)


```r
ir <- IRanges(start = c(2,6), width = 2)
aggregate(rl, ir, FUN = mean)
```

```
## [1] 1.0 2.5
```

The coverage() function counts, for each integer, how many ranges overlap the integer. Convert an IRanges to a Rle

```r
ir <- IRanges(start = 1:10, width = 3)
rl <- coverage(ir)
rl
```

```
## integer-Rle of length 12 with 5 runs
##   Lengths: 1 1 8 1 1
##   Values : 1 2 3 2 1
```



### RleList
An RleList is simply a list of Rle. Rle’s can also be constructed from GRanges.
This often involves RleList where each element of the list is a chromosome. 


```r
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
rl <- coverage(gr)
rl
```

```
## RleList of length 1
## $chr1
## integer-Rle of length 12 with 5 runs
##   Lengths: 1 1 8 1 1
##   Values : 1 2 3 2 1
```



## GenomicRanges - Lists


We will discuss GRangesList which is a list of GRanges

```r
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4, width = 3))
gL <- GRangesList(gr1 = gr1, gr2 = gr2)
gL
```

```
## GRangesList object of length 2:
## $gr1
## GRanges object with 4 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      *
##   [2]     chr1       2-4      *
##   [3]     chr1       3-5      *
##   [4]     chr1       4-6      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
## 
## $gr2
## GRanges object with 4 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr2       1-3      *
##   [2]     chr2       2-4      *
##   [3]     chr2       3-5      *
##   [4]     chr2       4-6      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

A number of standard GRanges functions work

```r
seqnames(gL)
```

```
## RleList of length 2
## $gr1
## factor-Rle of length 4 with 1 run
##   Lengths:    4
##   Values : chr1
## Levels(2): chr1 chr2
## 
## $gr2
## factor-Rle of length 4 with 1 run
##   Lengths:    4
##   Values : chr2
## Levels(2): chr1 chr2
```

```r
#Kind of xxapply
shift(gL, 10)
```

```
## GRangesList object of length 2:
## $gr1
## GRanges object with 4 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1     11-13      *
##   [2]     chr1     12-14      *
##   [3]     chr1     13-15      *
##   [4]     chr1     14-16      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
## 
## $gr2
## GRanges object with 4 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr2     11-13      *
##   [2]     chr2     12-14      *
##   [3]     chr2     13-15      *
##   [4]     chr2     14-16      *
##   -------
##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

findOverlaps works slightly different. For GRangesLists, we think of each element is a union of ranges. So we get an overlap if any range overlaps.

```r
findOverlaps(gL, gr2)
```

```
## Hits object with 4 hits and 0 metadata columns:
##       queryHits subjectHits
##       <integer>   <integer>
##   [1]         2           1
##   [2]         2           2
##   [3]         2           3
##   [4]         2           4
##   -------
##   queryLength: 2 / subjectLength: 4
```

Note how the queryLength is 2 and not 20. What we know from the first row of this output is that some range
in gL[[2]] overlaps the range gr[1].

This is actually a feature if we think of the GRangesList as a set of transcript, where each GRanges gives you the exon of the transcript. With this interpretation, 
findOverlaps tells you whether or not the transcript overlaps some region of interest, and this is true if any of the exons of the transcript overlaps the region.


# Biostrings





```r
BiocManager::install(c("Biostring"))
library(Biostrings)
```


The Biostrings package contains classes and functions for representing biological strings such as DNA,
RNA and amino acids. In addition the package has functionality for pattern matching (short read alignment) 
as well as a pairwise alignment function implementing Smith-Waterman local alignments and Needleman-Wunsch 
global alignments used in classic sequence alignment (see (Durbin et al. 1998)  for a description of these algorithms). 
There are also functions for reading and writing output such as FASTA files.


## Representing sequences

There are two basic types of containers for representing strings. One container represents a single string (say a chromosome or a single short read) and the other container represents a set of strings (say a set of short reads). 

There are different classes intended to represent different types of sequences such as DNA or RNA sequences.


```r
dna1 <- DNAString("ACGT-N")
dna1
```

```
##   6-letter "DNAString" instance
## seq: ACGT-N
```

```r
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2
```

```
##   A DNAStringSet instance of length 3
##     width seq
## [1]     4 ACGT
## [2]     4 GTCA
## [3]     4 GCTA
```
Note that the alphabet of a DNAString is an extended alphabet: - (for insertion) and N are allowed. 
In fact, IUPAC codes are allowed (these codes represent different characters, 
for example the code “M” represents either and “A” or a “C”). A list of IUPAC codes can be obtained by


```r
IUPAC_CODE_MAP
```

```
##      A      C      G      T      M      R      W      S      Y      K      V 
##    "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT"  "ACG" 
##      H      D      B      N 
##  "ACT"  "AGT"  "CGT" "ACGT"
```


Indexing into a DNAString retrieves a subsequence (similar to the standard R function substr), 
whereas indexing into a DNAStringSet gives you a subset of sequences.

```r
dna1[2:4]
```

```
##   3-letter "DNAString" instance
## seq: CGT
```

```r
dna2[2:3]
```

```
##   A DNAStringSet instance of length 2
##     width seq
## [1]     4 GTCA
## [2]     4 GCTA
```


Note that [[ allows you to get a single element of a DNAStringSet as a DNAString. 

```r
dna2[[2]]
```

```
##   4-letter "DNAString" instance
## seq: GTCA
```

DNAStringSet objects can have names, like ordinary vectors

```r
names(dna2) <- paste0("seq", 1:3)
dna2
```

```
##   A DNAStringSet instance of length 3
##     width seq                                               names               
## [1]     4 ACGT                                              seq1
## [2]     4 GTCA                                              seq2
## [3]     4 GCTA                                              seq3
```

The full set of string classes are

DNAString[Set]: DNA sequences.
RNAString[Set]: RNA sequences.
AAString[Set]: Amino Acids sequences (protein).
BString[Set]: “Big” sequences, using any kind of letter.

These classes seem very similar to standard characters() from base R, but there are important differences. 
The differences are mostly about efficiencies when you deal with either (a) many sequences or 
(b) very long strings (think whole chromosomes).

## Basic functionality

Basic character functionality is supported, like length, names.

c and rev (reverse the sequence).

width, nchar (number of characters in each sequence).

==, duplicated, unique.

as.charcater or toString: converts to a base character() vector.

sort, order.

chartr: convert some letters into other letters.

subseq, subseq<-, extractAt, replaceAt.

replaceLetterAt.

```r
width(dna2)
```

```
## [1] 4 4 4
```

```r
sort(dna2, decreasing=TRUE)
```

```
##   A DNAStringSet instance of length 3
##     width seq                                               names               
## [1]     4 GTCA                                              seq2
## [2]     4 GCTA                                              seq3
## [3]     4 ACGT                                              seq1
```

```r
rev(dna2)
```

```
##   A DNAStringSet instance of length 3
##     width seq                                               names               
## [1]     4 GCTA                                              seq3
## [2]     4 GTCA                                              seq2
## [3]     4 ACGT                                              seq1
```

```r
rev(dna1)
```

```
##   6-letter "DNAString" instance
## seq: N-TGCA
```


Note that rev on a DNAStringSet just reverse the order of the elements, whereas rev on a DNAString actually reverse the string.

## Biological functionality

There are also functions which are related to the biological interpretation of the sequences, including

reverse: reverse the sequence.

complement, reverseComplement: (reverse) complement the sequence.

translate: translate the DNA or RNA sequence into amino acids.


```r
translate(dna2)
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code,
## dna_codes[codon_alphabet], : in 'x[[1]]': last base was ignored
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code,
## dna_codes[codon_alphabet], : in 'x[[2]]': last base was ignored
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code,
## dna_codes[codon_alphabet], : in 'x[[3]]': last base was ignored
```

```
##   A AAStringSet instance of length 3
##     width seq                                               names               
## [1]     1 T                                                 seq1
## [2]     1 V                                                 seq2
## [3]     1 A                                                 seq3
```

```r
reverseComplement(dna1)
```

```
##   6-letter "DNAString" instance
## seq: N-ACGT
```

## Counting letters

We very often want to count sequences in various ways. Examples include:
  
Compute the GC content of a set of sequences.

Compute the frequencies of dinucleotides in a set of sequences.

Compute a position weight matrix from a set of aligned sequences.

There is a rich set of functions for doing this quickly.

alphabetFrequency, letterFrequency: Compute the frequency of all characters (alphabetFrequency) or only specific letters (letterFrequency).

dinucleotideFrequency, trinucleotideFrequency, 
oligonucleotideFrequeny: compute frequencies of dinucleotides (2 bases), trinucleotides (3 bases) and oligonucleotides (general number of bases).

consensusMatrix: consensus matrix; almost a position weight matrix.

Let’s look at some examples, note how the output expands to a matrix when you use the functions on a DNAStringSet:


```r
alphabetFrequency(dna1)
```

```
## A C G T M R W S Y K V H D B N - + . 
## 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0
```

```r
alphabetFrequency(dna2)
```

```
##      A C G T M R W S Y K V H D B N - + .
## [1,] 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [2,] 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## [3,] 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

```r
letterFrequency(dna2, "GC")
```

```
##      G|C
## [1,]   2
## [2,]   2
## [3,]   2
```

consensusMatrix() can be used to just compute the alphabet frequency for each position in the input sequences (most functions allows the return of probabilities with as.prob = TRUE):

```r
consensusMatrix(dna2, as.prob = TRUE)
```

```
##        [,1]      [,2]      [,3]      [,4]
## A 0.3333333 0.0000000 0.0000000 0.6666667
## C 0.0000000 0.6666667 0.3333333 0.0000000
## G 0.6666667 0.0000000 0.3333333 0.0000000
## T 0.0000000 0.3333333 0.3333333 0.3333333
## M 0.0000000 0.0000000 0.0000000 0.0000000
## R 0.0000000 0.0000000 0.0000000 0.0000000
## W 0.0000000 0.0000000 0.0000000 0.0000000
## S 0.0000000 0.0000000 0.0000000 0.0000000
## Y 0.0000000 0.0000000 0.0000000 0.0000000
## K 0.0000000 0.0000000 0.0000000 0.0000000
## V 0.0000000 0.0000000 0.0000000 0.0000000
## H 0.0000000 0.0000000 0.0000000 0.0000000
## D 0.0000000 0.0000000 0.0000000 0.0000000
## B 0.0000000 0.0000000 0.0000000 0.0000000
## N 0.0000000 0.0000000 0.0000000 0.0000000
## - 0.0000000 0.0000000 0.0000000 0.0000000
## + 0.0000000 0.0000000 0.0000000 0.0000000
## . 0.0000000 0.0000000 0.0000000 0.0000000
```










