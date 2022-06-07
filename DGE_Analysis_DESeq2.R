## install.packages("data.table")
## install.packages("BiocManager")
## install.packages("dplyr")
## install.packages("tibble")
## install.packages("DESeq2")
## BiocManager::install("ensembldb")
library(DESeq2)
library(ensembldb)

## Create Ensembl annotation object SQLite for reference gene using same reference GTF as featureCounts for consistency
## ensDbFromGtf("C:/ref/Mouse/Mus_musculus.GRCm38.100.gtf") ## only run of need to create sqlite file

ensdb<-EnsDb("C:/ref/Mouse/Mus_musculus.GRCm38.100.sqlite") ## Load sqlite for later use

d <- as.data.frame(data.table::fread("gene_count.txt")) ## data.table::fread load large file much faster
d <- d[,-2:-6]  ## Remove unnecessary column from gene count file for DESeq2 DGE
colnames(d) <- gsub(".bam", "", colnames(d))  ## Remove .bam from sample name
d <- d[rowMeans(d[-1] > 5,] ## Remove reads of gene that average to be below 5
d$symbol <- mapIds(ensdb,d$Geneid,"GENENAME","GENEID") ## Replace ensembl GENEID with more familiar gene symbol
d$Geneid <- NULL  ## Once the gene symbols are annotated the GENEID is no longer necessary (can be kepted if have other purpose)
d <- tibble::remove_rownames(d)  ## Remove rownames (if it exist to be replace with symbol)
d <- tibble::column_to_rownames(d,var = "symbol") ## This will place the symbol column into rownames

## This part will creates the coldata for the DGE
md <- data.frame(sample = colnames(d), dex = c(rep("exp",3),rep("ctrl",3))) ## This will create a data.frame with 2 columns: Sample ID, and comparison group.  The order can be change but for my purposes, the first 3 are experimentals and the last 3 are controls.

## Creating a dds object from everything we created
dds <- DESeqDataSetFromMatrix(countData = d,
                            colData = md,
                            design = ~dex,
                            tidy=F)
                          
vsd <- vst(dds) ## Creating a varianceStabilizingTransformation object which is used for PCA 
plotPCA(vsd,"dex") ## Examine sample variance of treatment and control group
 
dds <- DESeq(dds) ## Run default setting for differential expression analysis.  This function negates the need to run estimateSizeFactors, estimateDispersions, and nbinoWaldTest individually.

res <- as.data.frame(results(dds, contrast = c("dex","treatment","control")))  ## Create differential gene expression results from treatment_vs_control
