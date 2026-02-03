# Required inputs (not stored in git by default)

1) featurecounts_matrix.tsv
   Output of featureCounts (gene-level counts).

2) sampleinfo.tsv
   Tab-separated with columns:
   - sample
   - temperature
   - replicate (optional)

3) ros_targets.csv
   Two-column CSV:
   - GeneID
   - ROSTarget
   Where ROSTarget is one of: general, H2O2, O2-

4) (Optional) ros_gene_ids.txt
   A list of ROS-related gene IDs (e.g. 118 genes) to subset.
