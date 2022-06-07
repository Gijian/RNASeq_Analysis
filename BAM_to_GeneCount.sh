#!/bin/bash

## Simple script for extracting counts from all BAM files in ./bam/ folder
## Reference GTF file is downloaded from ensembl.org https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index (Newer version can be acquired) and stored on local drive
## Important to note that -s 2 is defined by the NEB Ultra II RNA Library Prep Kit which generate in a reverse strand format.  Other kits MAY not be the same.  0 is for unstranded, 1 is stranded, and 2 is reversed stranded.

featureCounts -T $(nproc) -p -s 2 -a "/mnt/c/ref/Mouse/Mus_musculus.GRCm38.100.gtf" -o genecount.txt ./bam/*.bam

## featureCounts
## Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. 
## Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PMID: 24227677.
