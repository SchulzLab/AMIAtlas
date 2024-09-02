This script generates the plots in the **supplementary figures 6 A**, using all cell types and with each cell type for **supplementary figures 6C - F**. The analysis is based on the SPONGE [1], to find the top 5 most connected genes, that have shared miRNAs targeting the genes. It uses the normalized miRNA expression and gene expression, and NCBI gene annotations for conversion of the gene symbols.

First, `make_interactionn_table.R` is used to create the matrix of known miRNA gene interactions. The resulting matrix variable `mat_intr` is used by `sponge_eachCelltype.R` to have the resultant ceRNA genes for each cell type.

Finally, to get the GO enriched functions for the ceRNA genes using all cell types for *molecular functions* and *biological processes* , we use `ceRNA_functional_barplots.R` that generates the bar plots in **supplementary figure 6 B.**

### Reference

1.  List, Markus, et al. "Large-scale inference of competing endogenous RNA networks with sparse partial correlation." Bioinformatics 35.14 (2019): i596-i604.
