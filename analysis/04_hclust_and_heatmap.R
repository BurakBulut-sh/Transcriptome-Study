#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript analysis/04_hclust_and_heatmap.R normalized_counts.tsv sampleinfo.tsv")
}

norm_file <- args[1]
sampleinfo_file <- args[2]

suppressPackageStartupMessages({
  library(tidyverse)
})

norm <- read.table(norm_file, header = TRUE, sep = "\t", check.names = FALSE)
meta <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

mat <- as.matrix(norm)
scaled <- scale(t(mat))   # scale per gene, samples as rows, then transpose back
scaled <- t(scaled)

# hclust (complete linkage) as in your notes
d <- dist(scaled)
hc <- hclust(d, method = "complete")

png("analysis/gene_hclust.png", width = 1200, height = 900)
plot(hc, labels = FALSE, main = "Gene hclust (complete)")
dev.off()

# Optional ComplexHeatmap if installed
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  ComplexHeatmap::png("analysis/heatmap_complexheatmap.png", width = 1200, height = 900)
  ComplexHeatmap::Heatmap(scaled, show_row_names = FALSE, show_column_names = TRUE)
  dev.off()
}

cat("Wrote: analysis/gene_hclust.png (and heatmap if ComplexHeatmap is available)\n")
