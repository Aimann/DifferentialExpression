## R script to run DESeq2 on a count matrix

# Load the command line arguments
args = commandArgs(trailingOnly=TRUE)
counts_file = args[1]
samples_file = args[2]
output_prefix = args[3]

# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.us.r-project.org")

BiocManager::install("DESeq2")
BiocManager::install("pheatmap")

# Load the required packages
library(DESeq2)
library(ggplot2)
library(pheatmap)

if (!dir.exists(output_prefix)){
  dir.create(output_prefix)
}

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
write
write.table(sizeFactors, file=paste0(output_prefix,"/","sizeFactors.txt"), sep="\t", quote=F, col.names='sizeFactor', row.names=T)

for (comb in comparisons){

    c1 = strsplit(comb, "_")[[1]][1]
    c2 = strsplit(comb, "_")[[1]][2]

    # Create a results table and adds the normalized counts to the results table
    res <- results(dds, contrast=c("condition", c2, c1), alpha=0.05, test='Wald', pAdjustMethod="BH")
    res <- cbind(res, as.data.frame(counts(dds, normalized=TRUE)))
    res <- res[order(res$padj), ]
    write.table(res, file=paste0(output_prefix,"/","deseq_results_", comb, ".txt"), sep="\t", quote=F, col.names=NA, row.names=T)

}

## Generates Correlation Plots & PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))
ggsave(paste0(output_prefix,"/","PCA.pdf"), width=10, height=10)

rld <- rlog(dds, blind=FALSE)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat, method="pearson")
heat <- pheatmap(rld_cor)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

save_pheatmap_pdf(heat, paste0(output_prefix,"/","Corr.pdf"), width=10, height=10)