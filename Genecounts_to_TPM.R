## The function used here was created by Kamil Slowikowski
## https://gist.github.com/slowkow/c6ab0348747f86e2748b

# install.packages("data.table")
# install.packages("BiocManager")
# BiocManager::install("ensembldb")
# BiocManager::install("AnnotationDbi")

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}
#######
## TPM function without meanFragmentLength requirement
counts_to_tpm_nofl <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

d <- as.data.frame(fread("genecount.txt"))

d <- d[,-2:-5]

colnames(d) <- gsub(".bam", "", colnames(d))

d <- tibble::column_to_rownames(tibble::remove_rownames(d), var = "Geneid"))

meanfraglength <- c(327,340,310,332,319,324) ## The mean fragment length was determined externally.  If this information is absence, use the previous function counts_to_tpm_nofl

tpm_count <- as.data.frame(counts_to_tpm(d[,-1], d[,1], meanfraglength)) ## Counts info is in the dataframe after first column.  Length information is in the first column. 

## tpm_count <- as.data.frame(counts_to_tpm_nofl(d[,-1], d[,1])) ## Use this if you do not have the information to mean fragment length

tpm_count$symbol <- AnnotationDbi::mapIds(ensembldb::EnsDb("C:/ref/Mouse/Mus_musculus.GRCm38.100.sqlite"),
                                  keys = row.names(tpm_count), 
                                  column = "GENENAME", 
                                  keytype = "GENEID") ## Annotation to gene symbol

write.csv(tpm_count, "TPM_genecounts.csv",row.names = F)

## AnnotationDbi
## Hervé Pagès, Marc Carlson, Seth Falcon and Nianhua Li (2021). AnnotationDbi: Manipulation of SQLite-based annotations in
## Bioconductor. R package version 1.54.1. https://bioconductor.org/packages/AnnotationDbi

## Other citations are in the "DGE_Analysis_DESeq2.R" 
