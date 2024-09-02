### Significance of the overlap of gene clusters and the miRNA target genes from miRNA clusters

The script `proportion_targets.R`, performs an overlap of significant genes (computed using maSigPro [1]) and the miRNA target gene clusters for each cell type, according to the workflow shown in **figure 4A**. Thereafter *Fisher’s exact test* is performed to find the significance of the overlap of the miRNA target genes and the genes in the cluster profiles. The code loops over different target gene correlation cut-offs' to find the best overlap.

> `Rscript proportion_targets.R`

The script generates dotplots as in figures 4B,C for CMs and FBs and Supplementary figures 4 A,B for ECs and HCs.

### Reference

1.  Nueda, María José, Sonia Tarazona, and Ana Conesa. "Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series." Bioinformatics 30.18 (2014): 2598-2602.
