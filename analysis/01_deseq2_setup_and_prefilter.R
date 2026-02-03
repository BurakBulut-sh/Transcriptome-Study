#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript analysis/01_deseq2_setup_and_prefilter.R featurecounts_matrix.tsv sampleinfo.tsv")
}

counts_file <- args[1]
sampleinfo_file <- args[2]

suppressPackageStartupMessages({
  library(DESeq2)
})

# featureCounts output includes annotation columns; keep only count columns
fc <- read.table(counts_file, header = TRUE, sep = "\t", quote = "", comment.char = "")
gene_id <- fc$Geneid
countmat <- fc[, -(1:6)]
rownames(countmat) <- gene_id

sampleinfo <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(sampleinfo) <- sampleinfo$sample

# Ensure column order matches sampleinfo
countmat <- countmat[, rownames(sampleinfo), drop = FALSE]

sampleinfo$condition <- factor(sampleinfo$temperature)
sampleinfo$condition <- relevel(sampleinfo$condition, ref = "15")

dds <- DESeqDataSetFromMatrix(countData = round(countmat), colData = sampleinfo, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

write.table(norm_counts, file = "analysis/data/normalized_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)
cat("Wrote: analysis/data/normalized_counts.tsv\n")
