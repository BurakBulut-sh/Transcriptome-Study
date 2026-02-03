# Workflow (RNA-seq → gene counts)

This workflow mirrors the project notes:
- Download SRA runs (prefetch → fasterq-dump)
- FastQC on paired FASTQ
- HISAT2 alignment to a pre-built index
- samtools view/sort
- featureCounts to generate gene-level counts

Configure everything in:
- workflow/config.env
- workflow/sra_run_ids.txt
