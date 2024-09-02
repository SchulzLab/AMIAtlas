### Generate miRNA expression dynamics as line plots


In order to check expression dynamics of known cell-type specific miRNAs compared between cell types (**figures** 3A,B,C), we use the script `generate_miRNA_expression_plots.R` as shown below:


>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-208a-5p -g miR-208a-5p -o mmu-miR-208a-5p`
>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-208a-3p -g miR-208a-3p -o mmu-miR-208a-3p`


>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-126b-3p -g miR-126b-3p -o mmu-miR-126b-3p`
>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-126b-5p -g miR-126b-5p -o mmu-miR-126b-5p`


>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-21a-3p -g miR-21a-3p -o mmu-miR-21a-3p`
>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-21a-5p -g miR-21a-5p -o mmu-miR-21a-5p`


### Generate miRNA expression heatmaps and functional enrichment of the expression clusters

The script `functional_enrichment_select_mirnas.R` generates the heatmaps of the normalized expression, for each cell type, generates clusters based on hierarchical clustering, and performs GO functional enrichment of the target genes for all miRNAs in cluster (**figures** 3D-G). The functions for each cell type and for each miRNA cluster are written in individual text files.
