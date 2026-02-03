# ROS × Temperature gradient in *Chironomus riparius*: RNA-seq + in vivo ROS (R + Slurm)

This repository contains:
1) A simple RNA-seq processing workflow (SRA download → FastQC → HISAT2 alignment → BAM sorting → featureCounts).
2) R analysis scripts for ROS-related gene expression (DESeq2 filtering/normalization, PCA, limma-based contrasts, volcano plots, hclust/heatmaps).

The scripts are based on the original project notes and were converted into reproducible, parameterized files.

## Associated paper
Bulut B. et al. (2025) *Ecology and Evolution* 15:e72625  
DOI: https://doi.org/10.1002/ece3.72625

## Data availability (from the paper)
- ENA accession: PRJEB89193
- Zenodo (primary data & metadata): https://doi.org/10.5281/zenodo.17184450

## What is NOT in this repo
- No raw FASTQ/BAM files are tracked.
- No large count matrices are tracked by default.

## Repository structure
- workflow/   : Slurm/Bash scripts to generate gene-level counts
- analysis/   : R scripts for downstream analysis and plots
- analysis/data/ : documentation of required input files

## Quick start (HPC)
1) Copy config:
   cp workflow/config.example.env workflow/config.env
   Edit workflow/config.env with your paths and module names.

2) Add run IDs:
   Edit workflow/sra_run_ids.txt

3) Run steps:
   bash  workflow/01_download_fastq.sh
   sbatch workflow/02_qc_fastqc.sbatch
   sbatch workflow/03_align_hisat2.sbatch
   bash  workflow/04_counts_featurecounts.sh

4) Analysis:
   Rscript analysis/01_deseq2_setup_and_prefilter.R analysis/data/featurecounts_matrix.tsv analysis/data/sampleinfo.tsv
   Rscript analysis/02_pca_antioxidant_genes.R      analysis/data/normalized_counts.tsv analysis/data/sampleinfo.tsv
   Rscript analysis/03_limma_volcano_plots.R        analysis/data/normalized_counts.tsv analysis/data/sampleinfo.tsv analysis/data/ros_targets.csv
   Rscript analysis/04_hclust_and_heatmap.R         analysis/data/normalized_counts.tsv analysis/data/sampleinfo.tsv
