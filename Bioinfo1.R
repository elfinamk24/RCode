## Differential Gene Expression Analysis with DESeq2

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")

library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load example count data (Replace with real dataset)
counts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("metadata.csv", row.names = 1)

# Ensure metadata matches count data
meta <- meta[match(colnames(counts), rownames(meta)), ]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)

# Pre-filtering to remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq analysis
dds <- DESeq(dds)
res <- results(dds)

# Order results by adjusted p-value
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), "deseq2_results.csv")

# Data visualization
# Volcano plot
res$logP <- -log10(res$padj)
ggplot(res, aes(x = log2FoldChange, y = logP)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differential Gene Expression")

# Heatmap of top 50 differentially expressed genes
topGenes <- rownames(head(res[order(res$padj), ], 50))
normCounts <- counts(dds, normalized = TRUE)
normCounts <- normCounts[topGenes, ]
pheatmap(normCounts, scale = "row", clustering_distance_rows = "correlation")
