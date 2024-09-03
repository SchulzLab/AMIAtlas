library(readr)
library(maSigPro)
library(MASS)
# ======================= get the normalized counts from DESEQ2 ===========================

normcounts_HC <- read.table(file.path("/projects/amitimeseries/work/salmon_run/data/normalized_genecounts_HC.txt"), header=TRUE)

# ==================== merge with the transcript-gene annotations ======================

tx2gene_mus = read_csv(file.path("/projects/amitimeseries/work/salmon_run/tx2gene.csv"), show_col_types = FALSE)
# normcounts_HC$txname = rownames(normcounts_HC)
# merged = merge(x = tx2gene_mus, y = normcounts_HC, all.y = TRUE, by.x = "TXNAME", by.y = "txname")


## ========================== build counts dataframe for maSigPro ==========================

# run for each cell-type
# data_ami = data.frame(); data_ami = merged[,c(-1)]; colnames(data_ami)[1] <- "gen"; rownames(data_ami) = merged[,1]#; View(data_ami)
# colnames(data_ami)

# ================= get the sample information ======================

edesign_HC <- read.table(file.path("samples_HC_diff.txt"), header=TRUE)# head(edesign_HC)
edesign_HC

dis_ami <- make.design.matrix(edesign_HC, degree = 5)#, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign_HC))) 


# ====================================== run maSigPro =======================================


# Performs a model comparison for each gene to detect genes with different trends in time course 
# experiments and applies maSigPro to the Isoforms belonging to selected genes.
# MyIso1 <- IsoModel(data=data_ami[,-1], gen=data_ami[,1], design=dis_ami, counts=TRUE, min.obs = 20)

# fit <- p.vector(data_ami[,-1], dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)
fit <- p.vector(normcounts_HC, dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)

#save the significant transcripts matrix
# write.table(fit$SELEC, file="sigprofiles_HC.txt", sep="\t", quote=F, col.names=NA)
tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)
get<-get.siggenes(tstep, vars="groups")

saveRDS(get, file = "MyIso_gene_HC.rds")

# readRDS(file = "MyIso_HC.rds")

## ======================= correlation =============
# my_data = read.table(file = "sigprofiles_HC.txt", header=TRUE, row.names = NULL, sep = "\t")

# colnames(my_data)


# row_sub = apply(my_data[, c(2:23)], 1, function(row) all(row !=0 ))
# my_data_noz = my_data[row_sub,]

# rownames(my_data_noz) = my_data_noz[, c(1)]   #View(my_data[, c(2)])

# my_data1 <- (my_data_noz[, c(2:23)])      #View(my_data[, c(2:23)])
# library(tidyverse)
# library("Hmisc")
# res2 <- rcorr(as.matrix(t(my_data1)))

# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# print("flattening ...")
# z = flattenCorrMatrix(res2$r, res2$P)
# yz3 = z[z$row %in% c("ENSMUST00000048129.6"), ]           #transcript id of piwi
# yz4 = z[z$column %in% c("ENSMUST00000048129.6"), ]           #transcript id of piwi
# print("beginning padjust ...")
# pvalue = c(yz3$p,yz4$p)
# z$padj = p.adjust(pvalue,method="BH")               #View(yz3)
# print("done this ...")
# yz = rbind(z[z$row %in% c("ENSMUST00000048129.6"), ], z[z$column %in% c("ENSMUST00000048129.6"), ])
# rowids = which(yz$column == c("ENSMUST00000048129.6"))
# yz[rowids,2] = yz[rowids,1]
# yz[rowids,1] = "ENSMUST00000048129.6"


# library(readr)
# tx2gene_mus = read_csv(file.path("../tx2gene.csv"), show_col_types = FALSE)
# merged_yz3 = merge(x = yz3, y = tx2gene_mus, by.x = "column", by.y = "TXNAME", all.x = TRUE); #View(merged_sigCM)

# # run biotype merging for each cell-type
# tx2biotype_mus = read_csv(file.path("../tx_biotype.csv"), show_col_types = FALSE)
# merged_diff_biotype_yz3 = merge(x = merged_yz3, y = tx2biotype_mus, by.x = "column", by.y = "tx", all.x = TRUE); 


# #run merging with gene symbol
# tx2symbol_mus = read_delim(file.path("../mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symb_yz3 = merge(x = merged_diff_biotype_yz3, y = tx2symbol_mus, by.x = "column", 
#                                  by.y = "Transcript stable ID version", all.x = TRUE);

# write.table(merged_diff_bio_symb_yz3, file="piwi_summary_HC.csv", sep="\t", quote=F, col.names=NA)


