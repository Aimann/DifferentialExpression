## R script to run DESeq2 on a count matrix

# Load the command line arguments
args = commandArgs(trailingOnly=TRUE)
counts_file = args[1]
samples_file = args[2]

# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.us.r-project.org")

BiocManager::install("DESeq2")

# Load the required packages
library(DESeq2)

# Set the working directory to the folder containing your count data
workdir <- getwd() 
setwd(workdir)

# Read count data from file
countData <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")

# Read sample metadata from file
metadata <- read.table(samples_file, header = TRUE, row.names = 1, sep = "\t")

# Create a DESeqDataSet object
dds <- DESeqDataSet(dds, colData = metadata, design= ~ assay + condition + assay:condition)

# Perform differential expression analysis
print("Performing LRT test using IP/Input design...")
dds <- DESeq(dds, test="LRT", reduced= ~ assay + condition)

## Extract the size factors
sizeFactors <- sizeFactors(dds)
write.table(sizeFactors, file="sizeFactors.txt", sep="\t", quote=F, col.names='sizeFactor', row.names=T)

## Adds the normalized counts to the results table
res <- results(dds)
res <- cbind(res, as.data.frame(counts(dds, normalized=TRUE)))
res <- res[order(res$padj), ]
write.table(res, file=paste0("deseq_results_LRT_enrichment.txt"), sep="\t", quote=F, col.names=NA, row.names=T)


