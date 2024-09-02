
### Cell type specific miRNAs through comparison of their targets


The script `cmoparing_targets.R` generates the lollipop plots for the top 5 cell type specific miRNAs (as in Figure 5B-E). This file requires the spearman correlations of the miRNAs with their targets. Also, to ascertain higher expression of a miRNA in the cell type of interest over others, the script requires the file containing miRNA normalized expressions.


To generate the functions of the cell type specific miRNAs (written in blue below each miRNA), we used `proportion_target.R`, that generates individual files for each miRNA, with their enriched functions.


