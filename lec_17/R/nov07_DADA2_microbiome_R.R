##***************************
## SOFTWARE TOOLS - DADA2 TUTORIAL
##
## Karl Cottenie
##
## 2024-11-06
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq"))

library(dada2)
library(phyloseq)

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

#### PART 1 - INTRODUCTION TO SCRIPT ----

#This file was updated from the notes of Sally Adamowicz

#INTRO: Today, we will go through a tutorial that highlights some of the important functions of the DADA2 package for 16S-based microbiome analysis. Please note that I did NOT write this script. During class, I will provide some guidance to help with understanding this script.

#This tutorial (with some additional commenting added in a few places) is from the following source. Please use this R file in conjunction with that source, as note I didn't copy over the text from the original source. Please cite the original source if you use or adapt this script.
#Tutorial by: Benjamin J Callahan, Joey McMurdie, Susan Holmes, October 24, 2023
#https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

#I also consulted and recommend:
#https://benjjneb.github.io/dada2/tutorial.html

#documentation
browseVignettes("dada2")

#### PART 2 - INTRODUCTION TO DADA2----

#Using example files that came with the package, system.file() will find the full names of files in packages
fnF1 <- system.file("extdata", "sam1F.fastq.gz", package = "dada2")
fnR1 <- system.file("extdata", "sam1R.fastq.gz", package = "dada2")

#How big is this? Small! The data are on hard disk. We can interact with zipped data stored on hard disk. So, R can be used for big data.
object.size(fnF1)

#By contrast, see:
x <- 1:100
object.size(x)
rm(x)

#setting up file names to house output data
filtF1 <- tempfile(fileext = ".fastq.gz")
filtR1 <- tempfile(fileext = ".fastq.gz")

#Note that the source data for this tutorial are Illumina paired-end sequencing data that are already demultiplexed (i.e. separate samples separated out using tags) and have had the primer sequences removed. Frequently, this step is provided by high-throughput sequencing facilities, but you should *ALWAYS CHECK* what steps have been performed on data you receive for analysis. DADA2 expects a separate fastq file per sample or two fastq files (one forward, and one reverse) per sample.

#Plotting quality of sequences, by sequence length. From the help for the plotQualityProfile function: "The plotted lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles."

#plot quality of forward sequences
plotQualityProfile(fnF1)

#plot quality of reverse sequences
plotQualityProfile(fnR1)

#What do you notice about the quality with cycle (i.e. with sequence length)? Drop off in quality. How fast is the dropoff in the forward set vs. the reverse set? e.g. What is the mean quality at a sequence length of 200 in the forward sequences vs. the reverse sequences?

#See documentation for this very helpful function"
?filterAndTrim

#filtering and trimming forward and reverse sequences. Trimming 10 nucleotides from the left (start) of all the sequences. Truncation of sequence length set to 240 for forward reads and 200 for reverse reads (i.e. trims off nucleotides after those sequence positions). Filtered out reads with any ambiguous nucleotides (by setting maxN = 0). Also, filtered out reads with more than two expected errors (by setting maxEE = 2).
filterAndTrim(fwd = fnF1, filt = filtF1, rev = fnR1, filt.rev = filtR1,
              trimLeft = 10, truncLen = c(240, 200), 
              maxN = 0, maxEE = 2,
              compress = TRUE, verbose = TRUE)

#Trimming and filtering criteria may need to be adjusted based on the quality of your run. As well, overlap between forward and reverse reads should be considered (which would depend upon the total length of the amplicon).

#Next comes "dereplication": collapsing sequences that are exactly the same.
derepF1 <- derepFastq(filtF1, verbose=TRUE)

derepR1 <- derepFastq(filtR1, verbose=TRUE)

class(derepF1)

derepF1

#As this is a small dataset for tutorial purposes, we can readily have a look at the data if we wish
derepF1$uniques
derepF1$quals
derepF1$map 

#Estimate error parameters for this dataset
errF <- learnErrors(derepF1, multithread = FALSE)
errR <- learnErrors(derepR1, multithread = FALSE)

#Infer sample composition
dadaF1 <- dada(derepF1, err = errF, multithread = FALSE)
dadaR1 <- dada(derepR1, err = errR, multithread = FALSE)
print(dadaF1)
print(dadaR1)

#Merging together forward and reverse reads, discarding reads that don't match
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose = TRUE)

#From the online tutorial: "The mergePairs(...) function returns a data.frame corresponding to each successfully merged unique sequence. The $forward and $reverse columns record which forward and reverse sequence contributed to that merged sequence."

#Remove chimeras. Chimeras are errors during PCR, whereby a read consists of part of one sequence and part of another. This can occur because of incomplete amplification, whereby the partial PCR product may prime a different sequence.
merger1.nochim <- removeBimeraDenovo(merger1, multithread = FALSE, verbose = TRUE)

class(merger1.nochim)

#Let's have a look
merger1.nochim

#Adding a second sample

# Assign filenames
fnF2 <- system.file("extdata", "sam2F.fastq.gz", package="dada2")
fnR2 <- system.file("extdata", "sam2R.fastq.gz", package="dada2")
filtF2 <- tempfile(fileext=".fastq.gz")
filtR2 <- tempfile(fileext=".fastq.gz")
# Filter and Trim
filterAndTrim(fwd=fnF2, filt=filtF2, rev=fnR2, filt.rev=filtR2, maxN=0, trimLeft=10, truncLen=c(240, 200), maxEE=2, compress=TRUE, verbose=TRUE)
## Read in 1500 paired-sequences, output 858 (57.2%) filtered paired-sequences.
# Dereplicate
derepF2 <- derepFastq(filtF2, verbose=TRUE)
## Dereplicating sequence entries in Fastq file: /tmp/RtmpZIeoAG/file20b92b76327d.fastq.gz
## Encountered 305 unique sequences from 858 total sequences read.
derepR2 <- derepFastq(filtR2, verbose=TRUE)
## Dereplicating sequence entries in Fastq file: /tmp/RtmpZIeoAG/file20b92df3419.fastq.gz
## Encountered 448 unique sequences from 858 total sequences read.
# Infer sample composition (using already learned error rates)
dadaF2 <- dada(derepF2, err=errF, multithread=FALSE)
## Sample 1 - 858 reads in 305 unique sequences.
dadaR2 <- dada(derepR2, err=errR, multithread=FALSE)
## Sample 1 - 858 reads in 448 unique sequences.
# Merge
merger2 <- mergePairs(dadaF2, derepF2, dadaR2, derepR2, verbose=TRUE)
## 731 paired-reads (in 6 unique pairings) successfully merged out of 852 (in 10 pairings) input.

#With that second sample processed, we can go ahead and create a sequence table.

#Finally we want to combine the inferred samples into one unified table. For that purpose we use makeSequenceTable:
seqtab <- makeSequenceTable(list(merger1, merger2))

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
## Identified 0 bimeras out of 7 input sequences.
dim(seqtab.nochim)
## [1] 2 7

#From the tutorial: "This is the final product of the dada2 pipeline, a matrix in which each row corresponds to a processed sample, and each column corresponds to an non-chimeric inferred sample sequence (a more precise analogue to the common "OTU table")." And, each cell contains the count of sequences for that variant. "From here we recommend proceeding forward with our friend the phyloseq package for further analysis."

#Additional comments. 

#Let's have a look in the viewer at object seqtab.nochim

#What does this look like? A lot like our "comm" object we have previously used for vegan! Now that we have the data in this format, we could use a variety of statistical methods to analyze the data. We could use the vegan package, for example. However, I would also recommend to consider using methods that consider the phylogenetic similarity among sequences, such as available in the phyloseq package.

#### PART 3 - CHALLENGE!----

#Use the ESV table from the above code and do something with it of interest to your group or on your own later.

#A few examples:

#Assignment taxonomy and generate biodiversity plots as per the following tutorial:
#https://benjjneb.github.io/dada2/tutorial.html

#Using functions available from the vegan package, build an individual based rarefaction curve for each sample. Are these communities well sampled? Or, would additional ESVs likely be discovered if sequencing depth were to be increased?

#Build a phylogenetic tree from the unique sequences and mark the relative abundances for the variants in each of the two samples.
#Example package for ideas the visualization: https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

