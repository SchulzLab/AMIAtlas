# BiocManager::install("SPONGE")

library(SPONGE)
library(ggplot2)
# =============== For the correct version of R =======================
# install.packages('devtools') 
# library(devtools)
# install_github('andreacirilloac/updateR')
# library(updateR)
# updateR()
## 

# ========== Prepare the miRNA expression table =================
setwd("~/spongeR_run/")
# sample_id x miMAT_ids

CM_miRNA_expression_alltimepoints <- read.delim("CM_miRNAExpression.txt", quote="")
FB_miRNA_expression_alltimepoints <- read.delim("FB_miRNAExpression.txt", quote="")
EC_miRNA_expression_alltimepoints <- read.delim("EC_miRNAExpression.txt", quote="")
HC_miRNA_expression_alltimepoints <- read.delim("HC_miRNAExpression.txt", quote="")

mirid_converted <- read.csv("./annotations/miRCarta - miRBase identifier conversion.csv", quote = "")

CM_anno = merge(x = CM_miRNA_expression_alltimepoints, y = mirid_converted, 
                by.x = "X",                # miRNA expression table
                by.y = "Query",            # conversion table
                all.y = TRUE)

EC_anno = merge(x = EC_miRNA_expression_alltimepoints, y = CM_anno, 
                by.x = "X",                # miRNA expression table
                by.y = "X",            # conversion table
                all.y = TRUE)

FB_anno = merge(x = FB_miRNA_expression_alltimepoints, y = EC_anno, 
                by.x = "X",                # miRNA expression table
                by.y = "X",            # conversion table
                all.y = TRUE)

HC_anno = merge(x = HC_miRNA_expression_alltimepoints, y = FB_anno, 
                by.x = "X",                # miRNA expression table
                by.y = "X",            # conversion table
                all.y = TRUE)

# drop the not required columns, assign the rownames to MIMAT ids 
rownames(HC_anno) = HC_anno$miRBase.accession

HC_anno_1 = subset(HC_anno, select = -c(X, miRBase.accession, d3_FB_2.mature,   
                                        d3_FB_3.mature, d28_HC_4 ))

#rename the columns

colnames(HC_anno_1) <- sub(".mature", "", colnames(HC_anno_1))
colnames(HC_anno_1) <- sub("ctrl_", "d0_", colnames(HC_anno_1))

all_celltypes_mirna_expr = HC_anno_1[complete.cases(HC_anno_1), ]
# ========== OBJECTIVE: sample_id x miMAT_ids
# convert to matrix and transpose 
all_anno_mat = as.matrix(t(all_celltypes_mirna_expr))

write.table(all_anno_mat, file="./expression/all_miRNAexpression_alltimepoints_anno.tsv", 
            sep = "\t", row.names = TRUE,
            quote = FALSE)

all_compiled_miRNA_expression_anno_1 <- read.delim("all_miRNAexpression_alltimepoints_anno.tsv", quote="", sep = "\t", row.names = 1)
all_compiled_miRNA_expression_anno = as.matrix(log(all_compiled_miRNA_expression_anno_1+10, base = 10))

View(all_compiled_miRNA_expression_anno)

# ========== Prepare the Gene expression table =================
#nrow(CM_gene_expression_anno)      16463
CM_gene_expression_anno <- read.delim("~/AMI_STEM_analysis/norm_counts_CM.csv", quote="", sep = ",")
#nrow(EC_gene_expression_anno)      18508
EC_gene_expression_anno <- read.delim("~/AMI_STEM_analysis/norm_counts_EC.csv", quote="", sep = ",")
#nrow(FB_gene_expression_anno)      18791
FB_gene_expression_anno <- read.delim("~/AMI_STEM_analysis/norm_counts_FB.csv", quote="", sep = ",")
#nrow(HC_gene_expression_anno)
HC_gene_expression_anno <- read.delim("~/AMI_STEM_analysis/norm_counts_HC.csv", quote="", sep = ",")

CM_EC_expression = merge(x = EC_gene_expression_anno, y = CM_gene_expression_anno, 
      by.x = "X",                # miRNA expression table
      by.y = "X",            # conversion table
      all.y = TRUE)

CM_EC_FB_expression = merge(x = FB_gene_expression_anno, y = CM_EC_expression, 
                         by.x = "X",                # miRNA expression table
                         by.y = "X",            # conversion table
                         all.y = TRUE)

CM_EC_FB_HC_expression = merge(x = HC_gene_expression_anno, y = CM_EC_FB_expression, 
                            by.x = "X",                # miRNA expression table
                            by.y = "X",            # conversion table
                            all.y = TRUE)

#save the NCBI ids file 
write.table(CM_EC_FB_HC_expression$X, file="all_NCBI_ids.tsv", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)

# get rid of rows with NAs
# all_celltypes_gene_expr = CM_EC_FB_HC_expression[complete.cases(CM_EC_FB_HC_expression), ]


#download the annotation file from ensembl---------- DONE

#import the downloaded annotation file, from DAVID

all_NCBI_anno <- read.delim("~/Downloads/conv_93C798CB7BF61622812016571.txt", quote="", sep = "\t")

# merge with gene_symbols

all_expression_anno = merge(x = all_NCBI_anno, y = CM_EC_FB_HC_expression, 
                               by.x = "From",
                               by.y = "X",            # conversion table
                               all.x = TRUE)

write.table(all_expression_anno, file="./significant_genes/compiled_anno.tsv", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)
# 
View(all_expression_anno)
nrow(all_expression_anno)
length(unique(all_expression_anno$To))


# Manually change the NCBI IDs that have same gene name

all_compiled_expression_anno <- read.delim("./significant_genes/compiled_anno.tsv", quote="", sep = "\t")
nrow(all_compiled_expression_anno)

# get rid of rows with NAs
all_celltypes_gene_expr = all_compiled_expression_anno[complete.cases(all_compiled_expression_anno), ]

rownames(all_celltypes_gene_expr) = all_celltypes_gene_expr$To   # check manually for duplicate gene names
colnames(all_celltypes_gene_expr)
all_celltypes_gene_expr_1 = subset(all_celltypes_gene_expr, select = -c(From, Species, Gene.Name,
                                                                        To))
# OBJECTIVE: sample_id x gene_symbols
colnames(all_celltypes_gene_expr_1)

all_gene_anno_mat = as.matrix(t(log(all_celltypes_gene_expr_1+10, base = 10)))
View(all_gene_anno_mat)


write.table(all_gene_anno_mat, file="./significant_genes/all_geneexpression_alltimepoints_anno.tsv", 
            sep = "\t", row.names = TRUE,
            quote = FALSE)
View(all_gene_anno_mat)
rownames(all_gene_anno_mat)
rownames(all_miRNAexpression_alltimepoints_anno)

identical(rownames(all_compiled_miRNA_expression_anno), rownames(all_gene_anno_mat))


# AMI: gene-miRNA interactions ==============
library(SPONGE)
#mat_intr obtained from make_interactionn_table.R
AMI_genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(gene_expr = all_gene_anno_mat, 
                                                                   mir_expr = all_compiled_miRNA_expression_anno, 
                                                                   mir_predicted_targets = mat_intr, #NULL,#targetscan_symbol
                                                                   log.level = "info",
                                                                   log.file = "sponge_run.log",
                                                                   # select.non.targets = TRUE,
                                                                   elastic.net = TRUE, 
                                                                   var.threshold = 0.5,
                                                                   coefficient.direction = "<",
                                                                   F.test = TRUE, 
                                                                   F.test.p.adj.threshold = 0.1,
                                                                   coefficient.threshold = -0.1
                                                                   )


View(AMI_genes_miRNA_candidates)

AMI_ceRNA_interactions <- sponge(gene_expr = all_gene_anno_mat,
                             mir_expr = all_compiled_miRNA_expression_anno,
                             mir_interactions = AMI_genes_miRNA_candidates)

View(AMI_ceRNA_interactions)

mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, 
                                            number_of_samples = nrow(all_gene_anno_mat))

sponge_plot_simulation_results(mscor_null_model)
AMI_ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = AMI_ceRNA_interactions, 
                                                   null_model = mscor_null_model)

View(AMI_ceRNA_interactions_sign)

AMI_ceRNA_interactions_fdr <- AMI_ceRNA_interactions_sign[which(AMI_ceRNA_interactions_sign$p.adj < 0.05),]
dim(AMI_ceRNA_interactions_fdr)
View((AMI_ceRNA_interactions_fdr))
install.packages("visNetwork")
sponge_plot_network(AMI_ceRNA_interactions_fdr, AMI_genes_miRNA_candidates)

network_centralities <- sponge_node_centralities(AMI_ceRNA_interactions_fdr)
ceRNA_interactions_fdr_weight <- AMI_ceRNA_interactions_fdr
ceRNA_interactions_fdr_weight$weight <- -log10(AMI_ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(AMI_ceRNA_interactions_fdr)
# install.packages("ggrepel")
sponge_plot_network_centralities(weighted_network_centralities, top = 10, base_size = 12)
sponge_plot_network_centralities(weighted_network_centralities, measure = "ev", top = 10, base_size = 12)
png("btw_plot.png", width = 784, height = 702)
sponge_plot_network_centralities(weighted_network_centralities, measure = "btw", top = 10, base_size = 12)
dev.off()
# Example: gene-miRNA interactions =======================================================
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mir_expr,
  mir_predicted_targets = targetscan_symbol)

genes_miRNA_candidates[1:2]
View(genes_miRNA_candidates)
# Example: ceRNA interactions ================ 
ceRNA_interactions <- sponge(gene_expr = gene_expr,
                             mir_expr = mir_expr,
                             mir_interactions = genes_miRNA_candidates)
head(ceRNA_interactions)
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, 
                                            number_of_samples = nrow(gene_expr))
a = sponge_plot_simulation_results(mscor_null_model)
ggsave("ceRNA_example_null_model.png", a)
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)

head(ceRNA_interactions_sign)
# ceRNA interaction network: per celltype ----------------------------------------------------
sponge_results_dir = "."
# for CM
load(paste0(sponge_results_dir, "ceRNAs/CM_SPONGE_network.RData"), ex1 <- new.env())
rm(list = ls(ex1) ) #ls(ex1) #ls.str(ex)
CM_AMI_ceRNA_interactions_sign = ex1$CM_AMI_ceRNA_interactions_sign
AMI_genes_miRNA_candidates = ex1$AMI_genes_miRNA_candidates

ceRNA_interactions_fdr <- CM_AMI_ceRNA_interactions_sign[
  which(CM_AMI_ceRNA_interactions_sign$p.adj < 0.5),]

sponge_plot_network(ceRNA_interactions_fdr, 
                    AMI_genes_miRNA_candidates)
network_analysis(ceRNA_interactions_fdr = ceRNA_interactions_fdr, celltype = "CM")
# for EC
load(paste0(sponge_results_dir, "ceRNAs/EC_SPONGE_network.RData"), ex <- new.env())  # load the object
ls(ex) #ls.str(ex)

EC_AMI_ceRNA_interactions_sign = ex$EC_AMI_ceRNA_interactions_sign
AMI_genes_miRNA_candidates = ex$AMI_genes_miRNA_candidates

ceRNA_interactions_fdr <- EC_AMI_ceRNA_interactions_sign[which(EC_AMI_ceRNA_interactions_sign$p.adj < 0.2),]

sponge_plot_network(ceRNA_interactions_fdr, 
                    AMI_genes_miRNA_candidates,
                    force.directed = TRUE)
network_analysis(ceRNA_interactions_fdr = ceRNA_interactions_fdr, celltype = "EC")

# for FB
rm(list = ls(ex) )
load(paste0(sponge_results_dir, "ceRNAs/FB_SPONGE_network.RData"), ex <- new.env())  # load the object
ls(ex) #ls.str(ex)

FB_AMI_ceRNA_interactions_sign = ex$FB_AMI_ceRNA_interactions_sign
AMI_genes_miRNA_candidates = ex$AMI_genes_miRNA_candidates

ceRNA_interactions_fdr <- FB_AMI_ceRNA_interactions_sign[which(FB_AMI_ceRNA_interactions_sign$p.adj < 0.5),]

sponge_plot_network(ceRNA_interactions_fdr, 
                    AMI_genes_miRNA_candidates,
                    force.directed = TRUE)
network_analysis(ceRNA_interactions_fdr = ceRNA_interactions_fdr, celltype = "FB")

# for HC
rm(list = ls(ex) )
load(paste0(sponge_results_dir, "ceRNAs/HC_SPONGE_network.RData"), ex <- new.env())  # load the object
ls(ex) #ls.str(ex)

HC_AMI_ceRNA_interactions_sign = ex$HC_AMI_ceRNA_interactions_sign
AMI_genes_miRNA_candidates = ex$AMI_genes_miRNA_candidates

ceRNA_interactions_fdr <- 
  HC_AMI_ceRNA_interactions_sign[which(HC_AMI_ceRNA_interactions_sign$p.adj < 0.5),]

sponge_plot_network(ceRNA_interactions_fdr, 
                    AMI_genes_miRNA_candidates,
                    force.directed = TRUE)
network_analysis(ceRNA_interactions_fdr = ceRNA_interactions_fdr, celltype = "HC")
# network analysis ----------------------------------------------------

network_analysis <- function(
  ceRNA_interactions_fdr,
  celltype
){
  network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
  
  ceRNA_interactions_fdr_weight <- ceRNA_interactions_fdr
  ceRNA_interactions_fdr_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
  weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
  
  cent = sponge_plot_network_centralities(weighted_network_centralities, top = 5)
  ggsave(paste0(sponge_results_dir, celltype, "_ceRNA_network_centralities.png"), cent)
  
  net_cent = sponge_plot_network_centralities(weighted_network_centralities, measure = "btw", top = 5)
  ggsave(paste0(sponge_results_dir, celltype, "_ceRNA_network_bw_centralities.png"), net_cent)

}





