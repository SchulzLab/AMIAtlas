### I
The script `gene_cluster_analysis.R` performs hierachical clustering on the normalized gene expression profiles on all four celltypes: cardiomyocytes (CMs), endothelial cells (ECs), fibroblasts (FBs) and hematopoietic cells (HCs).
The genes are chosen that are significantly different in expression at any timepoint (maSigpro results)
k was set to 10 clusters to highlight the expression dynamics of the significant genes across all timepoints (Figures 2:A,E,I,M) of the manuscript.


### II
The script `gene_cluster_analysis.R` also performs GO functional enrichments, comparing among the clusters that we identify from the previous step. It generates dot plots for select GO functional categories such as biological processes, or moleular functions.
, mentioned for each cell type.  
Some of these (Figures 2:B,F,J,N) of the manuscript.



