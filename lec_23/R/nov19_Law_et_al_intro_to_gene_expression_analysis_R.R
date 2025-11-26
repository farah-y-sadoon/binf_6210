######SOFTWARE TOOLS - CLASS 19 - INTRO TO GENE EXPRESSION ANALYSIS TUTORIAL BY LAW ET AL 2018----

#Authors: This tutorial on RNA-Seq analysis using Bioconductor packages is by: Charity Law, Monther Alhamdoosh, Shian Su, Xueyi Dong, Luyi Tian, Gordon K. Smyth and Matthew E. Ritchie. December 17, 2018. Some commenting added by Sally Adamowicz, last updated November 13, 2023.

#Please read through this tutorial at: https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

#The below script is from the above tutorial. Note that I did NOT copy in the commenting from the online tutorial here. So, you need to read over the tutorial online. Here, I added some supplementary commenting as well as some additional explorations of the data objects, so that we can feel more comfortable about understanding the data.

#I also recommend this set of slides about annotation resources:
#https://www.bioconductor.org/help/course-materials/2019/CSAMA/materials/lectures/lecture-08a-annotation.html#1

#The purpose here is to walk through a few features of the data and these packages, and then encourage you to go ahead and play with the data and also try out any Bioconductor vignettes that interest you!

####PACKAGES----

#If needed, install the below Bioconductor packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(c("limma", "Glimma", "edgeR", "Mus.musculus"))

#load packages
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

#If needed, install package R.utils and gplots from CRAN. Load packages.
#install.packages("R.utils")
library(R.utils)
#install.packages("gplots")
library(gplots)

####DATA----

#obtaining data. See online tutorial as well as the source paper (Sheridan et al. 2015).
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep="")) {
  R.utils::gunzip(i, overwrite=TRUE)}

#Let's have a look. Let's look at the first 5 rows of file #1.
read.delim(files[1], nrow=5)

#So, this line is combining the 9 files, taking only the first column (the gene IDs) and the third column (the counts) for each file. The data consist of three different types of cells, with three replicates each, hence 9 files. The function readDGE() will combine the 9 samples into a single object for analysis.
x <- readDGE(files, columns=c(1,3))

#what is the class?
class(x)

#Have a look at dimensions. There are 27179 genes and 9 samples. However, there are two elements in the DGEList object.
dim(x)
length(x)
names(x)

#So, our first element is a dataframe, and the second element is a matrix.
class(x[[1]])
dim(x[[1]])
class(x[[2]])
dim(x[[2]])

#TIP: I always like to have a look at a subset of the data, and I recommend that you do so too. It helps us to understand the contents of a particular data object.

#Let's have a look at the first element of our DGEList. This tells us the information about the samples and the total library size for each sample.
x[[1]]
class(x[[1]])

#This first data frame is called "samples". We can refer to the first element by name:
x$samples

#Now, let's have a look at the first fifty rows of the second element of our DGEList. This element is a dataframe which contains the gene IDs and counts of each gene for each sample. Note that we can use list-style indexing to have a look at the elements in this DGEList object.
x[[2]][1:100, ]

#Checking length of one column of this dataframe against our dimension information above for the DGEList object.
length(x[[2]][, 1])

#The second dataframe is called "counts". So, again we can refer to the second element of our DGEList by name if we wish, e.g.:
x$counts[1:10, ]

#Removing the GEO (Gene Expression Omnibus: https://www.ncbi.nlm.nih.gov/geo/) sample IDs, for simplicity for downstream analysis.

#First, let's look at colnames, so that we can clearly see what the below step is doing.
colnames(x)

#Authors of the tutorial wanted to clean up the sample names by omitting the GEO sample ids. This can be done by starting the sample names at character 12.
samplenames <- substring(text = colnames(x), first = 12, last = nchar(colnames(x)))
samplenames

#EXAMPLE: I am doing this a different way, as an example of using regular expressions to do this. This would be more flexible if there were variable numbers of characters.
testsamplenames <- sub(pattern = "^[A-Za-z]+[0-9]+_", x = colnames(x), replacement = "")

testsamplenames

all.equal(samplenames, testsamplenames)

#NOTE: The authors of this tutorial did this string processing by position. If we had variable numbers of characters at the start, this would be an example where we could look for the first underscore using a regular expression. Here, removing the starting information by position worked fine.

#setting up group information
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))

#creating a new variable for the groups (i.e. the three types of cells in this analysis)
x$samples$group <- group

#creating a new variable for the lane (on the Illumina sequencer) - it is good to keep information that could have influenced the experiment.
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

#retrieve annotations for mouse (see online tutorial for suggestions for other organisms).
geneid <- rownames(x)

#what do these gene id's look like? These are Entrez gene ids and can be used to retrieve other information of interest.
head(geneid)

#retrieving gene symbol and chromosome information, using the Entrez gene ids.
genes <- select(x = Mus.musculus, keys = geneid, columns = c("SYMBOL", "TXCHROM"), keytype = "ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

#adding gene annotations to our main data object
x$genes <- genes

#let's have a look
x

####DATA PREPROCESSING----

#Note that this analysis is looking for DIFFERENCES in gene expression among cell types for each gene, not absolute levels. So, note that the below metrics do not consider differences in gene length (which would impact the number of reads mapping to a particular gene).

#converting counts to counts per million and log counts per million
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#See tutorial for explanation. This is showing what the L parameter is for the formula used for the log method. Basically, the point of the formula is to avoid creating spurious appearance of large fold changes in gene expression for those genes with very low counts. Shrinks inter-sample log-fold changes towards zero.
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

summary(lcpm)

#removing genes lowly expressed

#first, summarizing unexpressed genes. 19% of genes in this dataset have a zero count across all 9 samples. (By contrast, we would wish to include genes that are expressed in one condition but not another.)
table(rowSums(x$counts==0)==9)

#filtering unexpressed and lowly expressed genes as well as genes expressed in few samples. The number of reads is considered as well as library size and group size.
keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#producing a figure to explore the effects of filtering
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#normalizing gene expression distributions, using TMM method
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#Changing the expression levels to provide an example of the impact of normalization
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#unsupervised clustering of samples
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

#This plot displays the samples in dimensions 1 and 2 (plot A) and 3 and 4 (plot B).

#interactive plot using Glimma package
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=TRUE)

####DIFFERENTIAL GENE EXPRESSION ANALYSIS----

#creating design matrix and contrasts
design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

#removing heteroscedasticity from count data
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#examining number of differentially expressed genes
summary(decideTests(efit))

#applying criterion of log-fold changes
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

head(tfit$genes$SYMBOL[de.common], n=20)

#Venn diagram of numbers of shared differentially expressed genes among conditions (three cell types)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

#writing output of differential gene expression results to file.
write.fit(tfit, dt, file="results.txt")

#examining individual gene results, starting with those with smallest (most significant) p-value
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)

head(basal.vs.ml)

#graphical representations of differential gene expression results
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

#interactive plot. Note I set launch to TRUE. You could set to FALSE of you don't want this to launch every time.
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)

basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

#see "camera method" at end of online tutorial

