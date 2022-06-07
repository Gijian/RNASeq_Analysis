## install.packages("data.table")
## install.packages("BiocManager")
## install.packages("dplyr")
## install.packages("tibble")
## install.packages("DESeq2")
## BiocManager::install("ensembldb")
library(DESeq2)

## Create Ensembl annotation object SQLite for reference gene using same reference GTF as featureCounts for consistency
# ensembldb::ensDbFromGtf("C:/ref/Mouse/Mus_musculus.GRCm38.100.gtf") ## only run of need to create sqlite file for gene symbol mapping later

ensdb<-ensembldb::EnsDb("C:/ref/Mouse/Mus_musculus.GRCm38.100.sqlite") ## Load sqlite

d <- as.data.frame(data.table::fread("gene_count.txt")) ## data.table::fread load large file much faster
d <- d[,-2:-6]  ## Remove unnecessary column from gene count file for DESeq2 DGE
colnames(d) <- gsub(".bam", "", colnames(d))  ## Remove .bam from sample name
d <- d[rowMeans(d[-1]) > 5,] ## Remove reads of gene that average to be below 5
d <- tibble::remove_rownames(d)  ## Remove rownames (if it exist to be replace with symbol)
d <- tibble::column_to_rownames(d,var = "Geneid") ## This will place the symbol column into rownames

## This part will creates the coldata for the DGE
md <- data.frame(sample = colnames(d), dex = c(rep("treatment",3),rep("control",3))) ## This will create a data.frame with 2 columns: Sample ID, and comparison group.  The order can be change but for my purposes, the first 3 are experimentals and the last 3 are controls.

## Creating a dds object from everything we created
dds <- DESeqDataSetFromMatrix(countData = d,
                            colData = md,
                            design = ~dex,
                            tidy=F)
                          
vsd <- vst(dds) ## Creating a varianceStabilizingTransformation object which is used for PCA 
plotPCA(vsd,"dex") ## Examine sample variance of treatment and control group
 
dds <- DESeq(dds) ## Run default setting for differential expression analysis.  This function negates the need to run estimateSizeFactors, estimateDispersions, and nbinoWaldTest individually.

res <- as.data.frame(results(dds, contrast = c("dex","treatment","control")))  ## Create differential gene expression results from treatment_vs_control
res$symbol <-  mapIds(x = ensembldb::EnsDb("C:/ref/Mouse/Mus_musculus.GRCm38.100.sqlite"), ## Mapping gene symbol from ensembl id using the same gtf files done from featureCounts
                      keys = row.names(res),
                      column = "GENENAME",
                      keytype = "GENEID") ## Replace ensembl GENEID with more familiar gene symbol

## Alternatively, I also used l2fcshrinkage to account for LFC of genes with low expression or high variation, which is common in sample size less than 5 per groups.
## I went with apeglm (Approximate Posterior Estimation for General Linear Model) approach to shrinkage which requires the installation of the "apeglm" package

# install.packages("apeglm")

res.ape <- as.data.frame(lfcShrink(dds = dds, coef = 2, type = "apeglm", apeMethod = "general"))
res.ape$symbol <- mapIds(x = ensembldb::EnsDb("C:/ref/Mouse/Mus_musculus.GRCm38.100.sqlite"),
                      keys = row.names(res),
                      column = "GENENAME",
                      keytype = "GENEID") 
## apeMethod = "general" is a bit slower but I prefer it.  The default can still be used.
## coef = 2 is used because it is the index designated for "treatment_vs_control" when you check with resultsNames(dds)

## Write out the DGE results as csv files for manual examination in excel examination and data visualization 
write.csv(res, "Results_DGE_of_Treatment_vs_Control.csv")
write.csv(res.ape, "Results_DGE_of_Treatment_vs_Control_shrinkage.csv")
                         
## DESeq2 normalized counts can be extracted out below for data visualization as well but I prefer TPM calculation.  
write.csv(as.data.frame(counts(dds,normalized = T)),"DESeq2_normalized_gene_counts.csv")

## tibble
## Kirill Müller and Hadley Wickham (2021). tibble: Simple Data Frames. R package version 3.1.6.
## https://CRAN.R-project.org/package=tibble

## data.table
## Matt Dowle and Arun Srinivasan (2021). data.table: Extension of `data.frame`. R package version 1.14.2.
## https://CRAN.R-project.org/package=data.table                         

## BiocManager
## Martin Morgan (2021). BiocManager: Access the Bioconductor Project Package Repository. R package version 1.30.16.
## https://CRAN.R-project.org/package=BiocManager
                         
## ensembldb
## Rainer J, Gatto L, Weichenberger CX (2019) ensembldb: an R package to create and use Ensembl-based annotation resources.
## Bioinformatics. doi:10.1093/bioinformatics/btz031

## DESeq2
## Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

## apeglm
## Zhu A, Ibrahim JG, Love MI (2018). “Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.” 
## Bioinformatics. doi: 10.1093/bioinformatics/bty895.
       
