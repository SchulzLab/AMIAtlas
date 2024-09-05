### Disease enrichment of orthologous cardiovascular disease genes

This script `run_disgenet.R` is based on the workflow given in figure 6A, which generates the dot-plot in figure 6B.

Enrichment is computed on the human orthologues of the differentially expressed genes for each cell type, and all mouse genes as the background. We use DOSE R library package to compute the enriched disease terms, esp. terms related to cardiovascular diseases. The enrichments are then compared between the cell types, to generate the resulting dot-plot. The file generated is stored in your results directory, at the same level as your other input data folders.

>Rscript run_disgenet.R


### Prerequisites


#### R Packages

`clusterprofiler`, `DOSE` for disease functional enrichments; `tidyr`, `dplyr` for data manupulations, `ggplot2` for plotting in R

#### Data

Reads the masigpro results with the significant genes for each celltype from the folder `significant_genes/<CT>_summary_genes.txt`, where <CT> is the cell type, CM, EC, FB or HC.