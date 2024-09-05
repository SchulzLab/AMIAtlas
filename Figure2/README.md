## Gene expression dynamics and their functional enrichment analysis

### Heatmap of gene expression profiles

The script `gene_cluster_analysis.R` performs hierachical clustering on the normalized gene expression profiles on all four celltypes: cardiomyocytes (CMs), endothelial cells (ECs), fibroblasts (FBs) and hematopoietic cells (HCs), and displays it as a heatmap. The pdf is stored in the `results` folder. The genes are chosen such that they are significantly different in their expression at any timepoint (using maSigpro). k was set to 10 clusters to highlight the expression dynamics of the significant genes across all timepoints (Figures 2:A,E,I,M) of the manuscript.

> Rscript gene_cluster_analysis.R


### Prerequisites


#### Package

`maSigPro` library to run significance analysis for the miRNAs; `ggplot2` for plotting, `pheatmap` library to generate heatmap and the clusters; `clusterprofiler` for functional annotations; `dplyr`, `tidyr` for dataframe manipulations in R.

#### Data

maSigPro generated significant genes are provided in the folder `significant_genes/<CT>_significant_genes_anno.txt` for CMs and ECs or `significant_genes/<CT>_significant_genes_remout_anno.txt` for the cell types with outlier samples i.e. FBs and HCs.


### Functional enrichment of gene clusters

The script `gene_cluster_analysis.R` also performs GO functional enrichment using the genes in each cluster, generated as a text file for each cell type.

The enrichment of the functions are then compared among the clusters identified from the previous step. It generates dot plots for *select* GO functions, mentioned for each cell type (**Figures** 2:B,F,J,N) of the manuscript for CMs, ECs, FBs and HCs.

> Rscript gene_cluster_analysis.R

### Gene expression plots 

In order to check expression dynamics of specific genes in a cell-type (**figures** 2C,D,G,H,O,P), we use the script `generate_tx_expression_plots.R` as shown below:

##### CM
>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000029580 -g Actb -o Actb`

>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000057329 -g Bcl2 -o Bcl2`

##### EC
>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000001847 -g Rac1 -o Rac1`

>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000006699 -g Cdc42 -o Cdc42`

##### FB
>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000026043  -g Col3a1 -o Col3a1` 

>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000062006 -g Rpl34 -o Rpl34`

##### HC
>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000027398 -g Il1b -o Il1b`

>`Rscript generate_tx_expression_plots.R -e ENSMUSG00000032725 -g Folr2 -o Folr2`

