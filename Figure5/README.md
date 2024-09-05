
### Cell type specific miRNAs through comparison of their targets


The script `cmoparing_targets.R` generates the lollipop plots for the top 5 cell type specific miRNAs (as in Figure 5B-E). This file requires the spearman correlations of the miRNAs with their targets. Also, to ascertain higher expression of a miRNA in the cell type of interest over others, the script requires the file containing miRNA normalized expressions.


To generate the functions enriched of the cell-type specific miRNAs (written in blue below each miRNA in the figure 5B, C, D, E), we used 

>`Rscript proportion_target.R`

that generates individual files for each miRNA, with their enriched functions (Supplementary table S3). 

### Prerequisites

#### Packages

Requires R packages: `gprofiler2` and `clusterProfiler`

#### Data

Data required for the script are provided in zenodo. The script requires miRNA spearman correlations, for each celltype: provided in the folder `correlations/<CT>_spearman_alltimepoints_anno.tsv`. CT stands for either CM, FB, EC or HC.
