### Generate miRNA expression dynamics as line plots


In order to check expression dynamics of known cell-type specific miRNAs compared between cell types (**figures** 3A,B,C), we use the script `generate_miRNA_expression_plots.R` as shown below:


### Prerequisites



#### Packages

`optparse` for parsing options; `tidyr` for data manipulation and `ggplot2` for plotting.


#### Data

miRNA expression files provided in the folder `expression/<CT>_mature_normalized_CPM.1.txt`. These are the CPM normalized miRNA expression, where CT stands for either CM, FB, EC or HC.




>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-208a-5p -g miR-208a-5p -o mmu-miR-208a-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-208a-3p -g miR-208a-3p -o mmu-miR-208a-3p`


>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-126b-3p -g miR-126b-3p -o mmu-miR-126b-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-126b-5p -g miR-126b-5p -o mmu-miR-126b-5p`


>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-21a-3p -g miR-21a-3p -o mmu-miR-21a-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-21a-5p -g miR-21a-5p -o mmu-miR-21a-5p`


### Generate miRNA expression heatmaps and functional enrichment of the expression clusters

The script `functional_enrichment_select_mirnas.R` generates the heatmaps of the normalized expression, for each cell type, generates clusters based on hierarchical clustering, and performs GO functional enrichment of the target genes for all miRNAs in cluster (**figures** 3D-G). The functions for each cell type and for each miRNA cluster are written in individual text files. The aggregated files are provided as **Supplementary table S1**. The arguments c referes to the celltype and k referes to the number of clusters. k mentioned in the following command refers to as used in the manuscript 

>`Rscript functional_enrichment_select_mirnas.R -c CM -k 7`

>`Rscript functional_enrichment_select_mirnas.R -c EC -k 5`

>`Rscript functional_enrichment_select_mirnas.R -c HC -k 6`

>`Rscript functional_enrichment_select_mirnas.R -c FB -k 6`


### Prerequisites



#### Packages

`ReactomePA`, `clusterProfiler` for functional enrichments; `maSigPro` for significance analysis; `pheatmap` for generating the heatmap and clustering.

#### Data

(i) miRNA spearman correlations, for each celltype: provided in the folder (in zenodo) `correlations/<CT>_spearman_alltimepoints_anno.tsv`. where CT stands for either CM, FB, EC or HC.

(ii) sample configurations as in the folders `configurations/`

(iii) miRNA expression files provided in the folder `expression/<CT>_mature_normalized_CPM.1.txt`. These are the CPM normalized miRNA expression, where CT stands for either CM, FB, EC or HC.

