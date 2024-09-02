### Disease enrichment of orthologous cardiovascular disease genes

This script is based on the workflow given in figure 6A, which generates the dot-plot in figure 6B.

Enrichment is computed on the human orthologues of the differentially expressed genes for each cell type, and all mouse genes as the background. We use DOSE R library package to compute the enriched disease terms, esp. terms related to cardiovascular diseases. The enrichments are then compared between the cell types, to generate the resulting dot-plot.
