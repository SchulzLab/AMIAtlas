setwd("~/spongeR_run/")

# conserved_interactions_mmu <- read.delim("~/spongeR_run/targetscan_interactions_mmu.txt")#larger set including non-conserved 
conserved_interactions_mmu <- read.delim("~/spongeR_run/conserved_interactions_mmu.txt")
dim(conserved_interactions_mmu)
#attach the MIMAT ids to these

write.table(unique(conserved_interactions_mmu$miRNA), file="conserved_mirnas_in_interaction.tsv", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)

miRNA_conversion = read.csv("~/Downloads/miRCarta - miRBase identifier conversion interaction.csv", quote="")
dim(miRNA_conversion)
# miRNA_conversion = read.csv("~/Downloads/miRCarta - miRBase identifier conversion.csv", quote="")#larger set


interactions_mmu = merge(x = miRNA_conversion, y = conserved_interactions_mmu, 
                               by.x = "Query",                # miRNA expression table
                               by.y = "miRNA",            # conversion table
                               all.x = TRUE)

View(interactions_mmu) # dim(interactions_mmu)
length(unique(interactions_mmu$miRBase.accession))
length(unique(interactions_mmu$Gene))

rm(mat_intr)
mat_intr = matrix(0:0, nrow = 10991, ncol = 370)
rownames(mat_intr) = unique(interactions_mmu$Gene)
colnames(mat_intr) = unique(interactions_mmu$miRBase.accession)
View(mat_intr)


for(k in 1:nrow(interactions_mmu)){
  mat_intr[interactions_mmu$Gene[k], interactions_mmu$miRBase.accession[k]] = 1
}
View(mat_intr)  
 
dim(mat_intr)  
# to check if the targetscan table (from sponge) has all miRNAs or not ======================
library(SPONGE)
View(mir_expr)
dim(mir_expr)
length(colnames(mir_expr))
class((targetscan_symbol) )
View(targetscan_symbol)     # 12238  x 348
