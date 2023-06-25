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
groups <- c(unique(metadata$condition))
comparisons = apply(combn(groups,2),2,paste,collapse='_')

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData, colData = metadata, design = ~ condition)

# Perform differential expression analysis
print("Performing Wald test...")
dds <- DESeq(dds)

## Extract the size factors
sizeFactors <- sizeFactors(dds)
write.table(sizeFactors, file="sizeFactors.txt", sep="\t", quote=F, col.names='sizeFactor', row.names=T)

for (comb in comparisons){

    c1 = strsplit(comb, "_")[[1]][1]
    c2 = strsplit(comb, "_")[[1]][2]

    # Create a results table and adds the normalized counts to the results table
    res <- results(dds, contrast=c("condition", c2, c1), alpha=0.05, test='Wald', pAdjustMethod="BH")
    res <- cbind(res, as.data.frame(counts(dds, normalized=TRUE)))
    res <- res[order(res$padj), ]
    write.table(res, file=paste0("deseq_results_", comb, ".txt"), sep="\t", quote=F, col.names=NA, row.names=T)
    
}


