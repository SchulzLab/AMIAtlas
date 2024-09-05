library(readr)
library(maSigPro)
library(MASS)
# ======================= get the normalized counts from DESEQ2 ===========================

normcounts_FB <- read.table(file.path("./expressions/normalized_genecounts_FB.txt"), header=TRUE)

# ================= get the sample information ======================

edesign_FB <- read.table(file.path("./configurations/samples_FB_diff.txt"), header=TRUE)# head(edesign_HC)

dis_ami <- make.design.matrix(edesign_FB, degree = 5)#, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign_FB))) 


# ====================================== run maSigPro =======================================


# Performs a model comparison for each gene to detect genes with different trends in time course 
# experiments and applies maSigPro to the Isoforms belonging to selected genes.
# MyIso1 <- IsoModel(data=data_ami[,-1], gen=data_ami[,1], design=dis_ami, counts=TRUE, min.obs = 20)

# fit <- p.vector(data_ami[,-1], dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)
fit <- p.vector(normcounts_FB, dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE)

#save the significant transcripts matrix
write.table(fit$SELEC, file="sig_gene_profiles_FB.txt", sep="\t", quote=F, col.names=NA)

tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)
get<-get.siggenes(tstep, vars="groups")

saveRDS(get, file = "MyIso_gene_FB.rds")

# readRDS(file = "MyIso_FB.rds")