#'   
#' ---
#' title: "case study 1: two group comparison"
#' author: "Veronique Storme"
#' date: "14/11/2018"
#' output: html_document
#' ---
#' 
#' 
#' bioconductor must have been installed as well as the edgeR package
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("BiocUpgrade")
#' biocLite("edgeR")
#' biocLite("limma")
#' 
#' The locfit package is also needed
#' biocLite("locfit")
#' 
#' other packages are needed for plotting and clustering, eg ggplo2
#' install.packages("ggplot2")
#' 
#' load the necessary libraries
#' 
#library(knitr)
library(limma)
library(edgeR)
library(locfit)
library(statmod)
library(Rcpp)
library(pheatmap)
library(gplots)
library(tidyr)

#' 
#' getting the version and correct citations
#' 
sessionInfo()
citation("base")
citation("edgeR")

#' 
#' this RNAseq experiment involves a control-treatment case
#' there are 3 independent biological samples for each treatment group
#' 
#' ### reading the counts from separate files
#' 
#' start by setting the working directory
#' 
## ------------------------------------------------------------------------

# working dir students:
setwd("//client/C$/MYCOURSES/RNAseq_v2/DEMO/data1")

# workin dir teacher:
setwd("//client/C$/Users/vesto/Documents/myCOURSES/RNAseq_v2/DEMO/data1")

#setwd("C:/MYCOURSES/RNAseq_v2/DEMO/data1")


#' 
#' counts are stored in 6 separate plain text files
#' 
## ------------------------------------------------------------------------
dir()

#' 
#' the file target_data1.txt gives the filename, the group and a brief description for each sample
#' 
## ------------------------------------------------------------------------
targets <- readTargets("target_data1.txt")
targets

#' 
#' read the tables of counts, calculate the sizes of the libraries and produce a DGElist object
#' 
## ------------------------------------------------------------------------
d <- readDGE(targets)

# the DGElist object contains 2 components: samples and counts
d$samples

# show the first 6 observations of the counts matrix
head(d$counts)

# number of genes and number of samples
dim(d$counts)

#' 
#' ### data exploration
#' 
#' (asinh is kind of log transformation but handles the zeroes better)
#' 
## ------------------------------------------------------------------------

# histogram of the counts after asinh transformation
hist(asinh(d$counts))

#' 
## ------------------------------------------------------------------------

# barplot of the library sizes: should be at least 20 million for arabidopsis
barplot(d$samples$lib.size*1e-6, names=1:6, ylab="Library size (millions)",cex.names=0.5)

#' 
#' get the counts in cpm (=counts divided by the libsize and multiplied by a million)
#' cpm(x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#' prior.count: average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE
#' ?cpm
#' 
## ------------------------------------------------------------------------

head(cpm(d))
?cpm

#' 
#' ### filtering
#' 
#' filter low expression tags (little power to detect DE for filtered tags)
#' 
#' a gene is required to have a count of 5-10 in a library to be considered expressed in that library
#' the expression must exceed 5-10 read counts in at 1 condition (if there are 3 replicates for each condition, then in at least 3 samples)
#' libsize divided by 1 million gives an idea about the relationship between raw count and cpm
#' thus 5 read counts correspond to 5 divided by the libsize and multiplied by a million
#' 
## ------------------------------------------------------------------------

# calculate nr of reads corresponding to 1 cpm = (library size/1000000)
prop <- d$samples$lib.size/1000000
prop 

# calculate cpm corresponding to 5 raw read counts
5/prop

#' 
#' this means that 5 read counts correspond to 0.34 to 1.12 cpm
#' filter on 1 cpm
#' 
## ------------------------------------------------------------------------

# keep the genes where the cpm is higher than 1 in at least 3 samples (assuming you have 3 observations/condition)
keep <- rowSums(cpm(d)>1) >=3
df <- d[keep,]

# number of genes after filtering
dim(df)

#' 
#' reset library sizes (must be done manually, using the colSums function)
#' 
## ------------------------------------------------------------------------
df$samples$lib.size <- colSums(df$counts)

#' 
## ------------------------------------------------------------------------

# verify that peak at zero counts has disappeared
hist(asinh(df$counts))

#' 
#' ### Normalization
#' 
#' apply TMM normalisation (TMM is the default method)
#' 
#' usage:
#' calcNormFactors(object, method=c("TMM","RLE","upperquartile","none"), refColumn = NULL,
#'       logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
#' refColumn: column (sample) to use as reference for method="TMM".
#' If refColumn is unspecified, the library whose upper quartile is closest to the mean upper quartile is used.
#' 
## ------------------------------------------------------------------------
df <- calcNormFactors(df)
df$samples

#' 
#' Note on cpm:
#' When TMM normalisation has been performed, then value for norm.factors, is then also used in cpm calulation
#' 
#' the effective library size = product of the original library
#' size (lib.size) and the scaling factor (norm.factors)
#' 
## ------------------------------------------------------------------------
lib.size.eff <- df$sample$lib.size*df$sample$norm.factors
lib.size.eff
# compare with
df$samples$lib.size
#' 

#' ### multidimensional scaling plot
#' 
#' The data can be explored by generating multi-dimensional scaling (MDS) plots.
#' This visualizes the differences between the expression profiles of different samples in two dimensions.
#' 
#' The distance between each pair of samples in the MDS plot is calculated as the leading fold change, defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples (default).
#' 
## ------------------------------------------------------------------------

?plotMDS
plotMDS(df)
plotMDS(df,labels = d$samples$description, cex=0.5)

#' 
#' ### dispersion estimation
#' 
#' > define the design of the experiment
#' 
## ------------------------------------------------------------------------

levels(df$samples$group)

# command in case you want to relevel the group levels
df$samples$group <- relevel(df$samples$group, ref="control")

design = model.matrix(~ group, data=df$samples)

design

# add rownames
rownames(design) <- df$samples$group
design
#' 
#' > non-robust estimate of the prior degrees of freedom
#' 
#' common, trended and tagwise dispersions are estimated
#' 
#' When robust=FALSE, there is only one prior.df estimated, when robust=TRUE, each gene gets a prior.df value 
#' Ref: Phipson et al. 2013 "empirical Bayes in the presence of exceptional cases with application to microarray data"
#' 
## ------------------------------------------------------------------------
df <- estimateDisp(df,design=design, robust=FALSE)
df$prior.df

#' 
#' > robust estimate of the prior degrees of freedom
#' 
#' robust=TRUE : highly variable genes are less squeezed to the trend, resulting in a lower estimate of prior.df (=weight given to the trend)
#' 
## ------------------------------------------------------------------------
dfr <- estimateDisp(df,design=design, robust=TRUE)
dfr$prior.df[1:5]

#' 
#' The following plot shows the biological coefficient of variation for each gene against its abundance with the estimated common dispersion and trended dispersions
#' 
## ------------------------------------------------------------------------
plotBCV(dfr,cex=0.4)

#' outliers are marked by low prior.df values (<1)
## ------------------------------------------------------------------------
outliers = dfr$counts[dfr$prior.df < 1,]
dim(outliers)
head(outliers)

# if there are outliers, save them as a csv file
write.csv(outliers, file="outliers.csv")

# A gene sign DE but present in the outliers list is most probably a false positive

#' 
#' ### fitting a quasi-likelihood model
#' 
## ------------------------------------------------------------------------
fit = glmQLFit(dfr, design, robust=TRUE)

#' 
#' ### Testing for DE genes
#' 
#' The model fitted is:
#' 
#' log(E(Y/efflib)) = b0 + b1X (X=1 when group=treated and 0 when group=control)
#' 
#' verify with
#' 
## ------------------------------------------------------------------------
design

#' 
#' H0: log(E(Y/efflib|treated) = log(E(Y/efflib|control) 
#' 
#' or log(E(Y/efflib|treated) - log(E(Y/efflib|control) = 0
#' 
#' log(E(Y/efflib)|treated) = b0 + b1
#' 
#' log(E(Y/efflib)|control) = b0
#' 
#' Thus log(E(Y/efflib|treated) - log(E(Y/efflib|control) = b1 = 0
#' 
#' Thus we need to test whether b1=0
#' 
## ------------------------------------------------------------------------

DE = glmQLFTest(fit, coef=2)

# this object is a list with a component table
head(DE$table)

#' Note: logFC is here log2FC
#' 
#' topTags() produces a table of the top differentially expressed tags (default n=10)
#' default adjustment for multiple testing with the BH method (FDR column)
#' sorted by (unadjusted) pvalue (default) or by absolute log2-fold change or none
#' 
## ------------------------------------------------------------------------

DE.fdr = topTags(DE,n=Inf,sort.by="none")
head(DE.fdr$table)

# number of DE genes
sum(DE.fdr$table$FDR<0.05)

#' export as csv file
#' 
## ------------------------------------------------------------------------
write.csv(DE.fdr$table,file="de_fdr_logFC.csv")

#' 
#' > use of decideTestsDGE function
#' 
#' total nr of DE genes at FDR < 0.05 is :
#' 
## ------------------------------------------------------------------------
de <- decideTestsDGE(DE, p=0.05, adjust="BH")
de[1:5,]
summary(de)

#' 
#' entries for -1, 0, and 1  are for down-regulated, nonDE and up-regulated tags, respectively
#' 
#' Note that with the previous non-robust nonQL method, 15 genes were found to be DE and upregulated
#' 
#' The QL method with the robust option has a better control of the FDR
#' 
#' ### using a different parameterization: the no intercept model
#' 
## ------------------------------------------------------------------------

group = c(rep("control",3),rep("treated",3))
design2 = model.matrix(~ 0 + group)
rownames(design2) <- group
design2

#' 
#' re-estimate dispersions
#' 
## ------------------------------------------------------------------------
dfr <- estimateDisp(df,design=design2, robust=TRUE)

#' 
#' re-fit the model
#' 
## ------------------------------------------------------------------------
fit2 = glmQLFit(dfr, design, robust=TRUE)

#' 
#' Now we have:
#' 
#' log(E(Y/efflib)) = b1X1 + b2X2 
#' (X1=1 when group=control, X2=1 when group=treated)

#' log(E(Y/efflib|control)) = b1
#' log(E(Y/efflib|treatment)) = b2
#' 
#' testing for DE
#' 
#' now we need a different contrast
#' 
## ------------------------------------------------------------------------
con <- makeContrasts(grouptreated - groupcontrol, levels=design2)
con
DE = glmQLFTest(fit, contrast=con)
head(DE$table)

DE.fdr = topTags(DE,n=Inf,sort.by="none")
head(DE.fdr$table)
