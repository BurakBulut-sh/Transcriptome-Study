#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript analysis/03_limma_volcano_plots.R normalized_counts.tsv sampleinfo.tsv ros_targets.csv")
}

norm_file <- args[1]
sampleinfo_file <- args[2]
targets_file <- args[3]

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

norm <- read.table(norm_file, header = TRUE, sep = "\t", check.names = FALSE)
meta <- read.table(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
targets <- read.csv(targets_file, header = TRUE, stringsAsFactors = FALSE)

# log2 transform with pseudocount 1 (as in your notes)
expr <- log2(as.matrix(norm) + 1)

# design: compare each temp vs 15
meta$temperature <- factor(meta$temperature)
meta$temperature <- relevel(meta$temperature, ref = "15")
design <- model.matrix(~ 0 + temperature, data = meta)
colnames(design) <- gsub("temperature", "", colnames(design))

fit <- lmFit(expr, design)

# contrasts vs 15°C
temps <- setdiff(levels(meta$temperature), "15")
contr <- makeContrasts(contrasts = paste0(temps, "-15"), levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

plot_one <- function(tt, title) {
  tt$GeneID <- rownames(tt)
  tt <- merge(tt, targets, by.x = "GeneID", by.y = "GeneID", all.x = TRUE)
  tt$ROSTarget[is.na(tt$ROSTarget)] <- "unknown"

  tt$neglog10FDR <- -log10(tt$adj.P.Val + 1e-300)

  ggplot(tt, aes(x = logFC, y = neglog10FDR, color = ROSTarget)) +
    geom_point(alpha = 0.8, size = 1.6) +
    theme_minimal() +
    labs(title = title, x = "log2 fold change", y = "-log10(FDR)")
}

dir.create("analysis/volcano", showWarnings = FALSE)

for (t in temps) {
  tt <- topTable(fit2, coef = paste0(t, "-15"), number = Inf, sort.by = "none")
  p <- plot_one(tt, paste0(t, "°C vs 15°C"))
  ggsave(filename = file.path("analysis/volcano", paste0("volcano_", t, "_vs_15.png")),
         plot = p, width = 7, height = 5, dpi = 300)
}

cat("Wrote volcano plots to analysis/volcano/\n")
