## R script to run DESeq2 on a count matrix

# Load the command line arguments
args = commandArgs(trailingOnly=TRUE)
deseq_file = args[1]
samples_file = args[2]
output_file = args[3]

# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.us.r-project.org")

BiocManager::install("DEGreport")

# Load the required packages
library(DEGreport)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
# Set the working directory to the folder containing your count data
workdir <- getwd() 
setwd(workdir)

# Read count data from file
deseq_results <- read.table(deseq_file, header = TRUE, row.names = 1, sep = "\t")

deseq_results <- drop_na(deseq_results)

## Filters for significant genes
deseq_results <- deseq_results[deseq_results$padj <= 0.05,]
# deseq_results <- head(deseq_results, 500)

deseq_results <- deseq_results[,7:ncol(deseq_results)]

## log2 normalize
deseq_results <- log2(deseq_results + 1)

metadata <- read.table(samples_file, header = TRUE, row.names = 1, sep = "\t")

res <- degPatterns(deseq_results, metadata, time="condition", col="genotype", minc=1)

resout <- res$df
write.table(resout, file=output_file, sep="\t", quote=F, col.names=NA, row.names=T)

degPlotCluster(res$normalized, "condition", "genotype")
ggsave("degPlotCluster.pdf", width=10, height=10)


