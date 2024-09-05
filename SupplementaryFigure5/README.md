

### Generate miRNA expression dynamic plots


In order to verify higher expression visually in one cell type of interest, compared to others (figures in **Supplementary figures 5A-D**), we use `generate_miRNA_expression_plots.R` as shown below, to generate the miRNA expression plots for all cell types:

### Prerequisites

#### Packages

`optparse` for parsing options; `tidyr` for data manipulation and `ggplot2` for plotting.


#### Data

miRNA expression files provided in the folder `./expression/<CT>_mature_normalized_CPM.1.txt`. These are the CPM normalized miRNA expression, where CT stands for either CM, FB, EC or HC.

##### CMs 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-465b-5p -g miR-465b-5p -o mmu-miR-465b-5p` 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-465c-5p -g miR-465c-5p -o mmu-miR-465c-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-3057-5p -g miR-3057-5p -o mmu-miR-3057-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-7068-3p -g miR-7068-3p -o mmu-miR-7068-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-3061-5p -g miR-3061-5p -o mmu-miR-3061-5p`

##### ECs

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-3470b-3p -g miR-3470b-3p -o mmu-miR-3470b-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-706 -g miR-706 -o mmu-miR-706` 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-466i-5p -g miR-466i-5p -o mmu-miR-466i-5p` 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-1291 -g miR-1291 -o mmu-miR-1291` 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-3470a -g miR-3470a -o mmu-miR-3470a`

##### HCs

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-466d-3p -g miR-466d-3p -o mmu-miR-466d-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-669f-3p -g miR-669f-3p -o mmu-miR-669f-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-297b-5p -g miR-297b-5p -o mmu-miR-297b-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-744-3p -g miR-744-3p -o mmu-miR-744-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-196a-5p -g miR-196a-5p -o mmu-miR-196a-5p`


##### FBs 

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-344d-3p -g miR-344d-3p -o mmu-miR-344d-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-129-5p -g miR-129-5p -o mmu-miR-129-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-34c-3p -g miR-34c-3p -o mmu-miR-34c-3p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-181b-5p -g miR-181b-5p -o mmu-miR-181b-5p`

>`Rscript generate_miRNA_expression_plots.R -e mmu-miR-493-3p -g miR-493-3p -o mmu-miR-493-3p`

