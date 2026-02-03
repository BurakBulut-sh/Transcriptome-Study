#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript analysis/02_pca_antioxidant_genes.R normalized_counts.tsv sampleinfo.tsv")
}

norm_file <- args[1]
sampleinfo_file <- args[2]

suppressPackageStartupMessages({
  library(ggplot2)
})

norm <- read.table(norm_file, header = TRUE, sep = "\t", check.names = FALSE)
sampleinfo <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# norm: genes x samples
mat <- as.matrix(norm)
mat_t <- t(mat)                 # samples x genes
pca <- prcomp(mat_t, scale. = TRUE)

scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)

meta <- sampleinfo
scores <- merge(scores, meta, by = "sample")

p <- ggplot(scores, aes(x = PC1, y = PC2, color = as.factor(temperature))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of RNA-seq data", color = "Temperature")

ggsave("analysis/pca_PC1_PC2.png", p, width = 7, height = 5, dpi = 300)
cat("Wrote: analysis/pca_PC1_PC2.png\n")
