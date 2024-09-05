library(readr)
library(maSigPro)
library(MASS)
# ============== get the normalized counts from DESEQ2 ==================

normcounts_CM <- read.table(file.path(
    "./expression/normalized_genecounts_CM.tsv"), header=TRUE)


# ============== merge with the transcript-gene annotations ================

# tx2gene_mus = read_csv(file.path("/work/salmon_run/tx2gene.csv"), show_col_types = FALSE)
# normcounts_CM$txname = rownames(normcounts_CM)
# merged = merge(x = tx2gene_mus, y = normcounts_CM, all.y = TRUE, by.x = "TXNAME", by.y = "txname")


## ========================== build counts dataframe for maSigPro ==========================

# run for each cell-type
# data_ami = data.frame(); data_ami = merged[,c(-1)]; colnames(data_ami)[1] <- "gen"; rownames(data_ami) = merged[,1]#; View(data_ami)
# colnames(data_ami)

# ================= get the sample information ======================

edesign_CM <- read.table(file.path("./configurations/samples_CM_diff.txt"), header=TRUE)# head(edesign_CM)
dis_ami <- make.design.matrix(edesign_CM, degree = 5)#, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign_CM)))      # rownames(edesign_ami);# colnames(data_ami)

# ========================= run maSigPro ========================


# Performs a model comparison for each gene to detect genes with different trends in time course 
# experiments and applies maSigPro to the Isoforms belonging to selected genes.
# MyIso1 <- IsoModel(data=data_ami[,-1], gen=data_ami[,1], design=dis_ami, counts=TRUE, min.obs = 20)

# fit <- p.vector(data_ami[,-1], dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)
fit <- p.vector(normcounts_CM, dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)

#save the significant transcripts matrix
write.table(fit$SELEC, file="sig_gene_profiles_CM.txt", sep="\t", quote=F, col.names=NA)

tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)
get<-get.siggenes(tstep, vars="groups")

saveRDS(get, file = "MyIso_gene_CM.rds")

# readRDS(file = "MyIso_CM.rds")




