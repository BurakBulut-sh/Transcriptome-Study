# R analysis scripts

These scripts implement the key analysis blocks from the project notes:
- Read featureCounts output
- Create sample info (temperature as factor; 15°C as reference)
- Pre-filter lowly expressed genes (>=10 counts in >=4 samples)
- PCA via prcomp (scaled, samples as rows)
- limma contrasts vs 15°C and volcano plots
- hclust of scaled expression + optional ComplexHeatmap

Inputs expected in analysis/data/ (see analysis/data/README.md).
