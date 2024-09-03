## ============================= on HC and EC long RNA seq data ==============================
setwd("/Users/ranjan/time_course_analysis/masigpro/")
library(readr)
library(maSigPro)
library(MASS)

library("optparse")

option_list = list(
    make_option(c("-c", "--celltype"), type="character", default=NULL, 
              help="cell type of interest, either CM, EC, FB or HC", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$c)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as the cell type", call.=FALSE)
} 


# =================== get the normalized counts from DESEQ2 =============================
rm(normcounts, normcounts_HC, normcounts_EC, normcounts_FB,normcounts_CM)
# normcounts_HC <- read.table(file.path("normalized_counts_HC.txt"), header=TRUE)
normcounts_HC <- read.table(file.path("normalized_genecounts_HC.txt"), header=TRUE)

# normcounts_EC <- read.table(file.path("normalized_counts_EC.txt"), header=TRUE)
normcounts_EC <- read.table(file.path("normalized_genecounts_EC.txt"), header=TRUE)

# normcounts_FB <- read.table(file.path("normalized_counts_FB.txt"), header=TRUE)
normcounts_FB <- read.table(file.path("normalized_genecounts_FB.txt"), header=TRUE)

# normcounts_CM <- read.table(file.path("normalized_counts_CM.txt"), header=TRUE)
normcounts_CM <- read.table(file.path("normalized_genecounts_CM.txt"), header=TRUE)

# normcounts = cbind(normcounts_HC, normcounts_EC, normcounts_FB, normcounts_CM)

# # =================== merge with the transcript-gene annotations ==========================
# tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)

# # run for each cell-type
# normcounts_HC$txname = rownames(normcounts_HC)				#colnames(normcounts)
# merged = merge(x = tx2gene_mus, y = normcounts_HC, all.y = TRUE, by.x = "TXNAME", by.y = "txname"); #View(normcounts)
# colnames(merged)
# ## ============== build counts dataframe for maSigPro (also for indv model) ==========================

# # run for each cell-type
# data_ami = data.frame(); data_ami = merged[,c(-1)]; colnames(data_ami)[1] <- "gen"; rownames(data_ami) = merged[,1]#; View(data_ami)
# colnames(data_ami)
# ## ===========  get the sample information and create the design matrix =====================
# edesign_HC <- read.table(file.path("samples_HC_diff.txt"), header=TRUE)        # head(edesign_HC); 

# edesign_EC <- read.table(file.path("samples_EC.txt"), header=TRUE)
# edesign_CM <- read.table(file.path("samples_CM.txt"), header=TRUE)
# edesign_FB <- read.table(file.path("samples_FB.txt"), header=TRUE)

# edesign_ami <- read.table(file.path("samples_full.txt"), header=TRUE)


# dis_ami <- make.design.matrix(edesign_CM, degree = 3)      # rownames(edesign_ami);# colnames(data_ami)
# dis_ami <- make.design.matrix(edesign_EC, degree = 4)      # rownames(edesign_ami);# colnames(data_ami)
# dis_ami <- make.design.matrix(edesign_FB, degree = 4)      # rownames(edesign_ami);# colnames(data_ami)
# dis_ami <- make.design.matrix(edesign_HC, degree = 5)#, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign_HC)))   
# rm(Myget1,MyIso1,see1)

# colnames(data_ami)
# ## =============================== run for IsoModel in masigpro ==============================


# # Performs a model comparison for each gene to detect genes with different trends in time course 
# # experiments and applies maSigPro to the Isoforms belonging to selected genes.
# MyIso1 <- IsoModel(data=data_ami[,-1], gen=data_ami[,1], design=dis_ami, counts=TRUE, min.obs = 6)

# # MyIso1 <- IsoModel(data=data_ami[1:20000,-1], gen=data_ami[1:20000,1], design=dis_ami, counts=TRUE, min.obs = 20)
# View(MyIso1$pvector$SELEC)     ; dim(MyIso1$pvector$SELEC)       #10748  matrix containing the expression values for significant genes
# View(MyIso1$pvector$p.vector)  ; dim(MyIso1$pvector$p.vector)       #vector containing the computed p-values, 18990
# View(MyIso1$pvector$p.adjusted)       # vector containing the adjusted p-values, 18990
# View(MyIso1$pvector$i)                #10748: number of significant genes
# View(MyIso1$pvector$g)  #18990; number of genes taken in the regression fit
# View(MyIso1$pvector$G)  #25980


# x_df  = as.data.frame(MyIso1$pvector$SELEC); x_df$tx = rownames(MyIso1$pvector$SELEC)
# y_df = as.data.frame(cbind(MyIso1$pvector$p.vector, MyIso1$pvector$p.adjusted)); y_df$tx = rownames(y_df)

# sig_CM = merge(x = x_df, y = y_df, by.x = "tx", by.y="tx", all.x=TRUE)

# View(sig_CM)

# # run geneID merging for each cell-type
# tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# merged_sigCM = merge(x = sig_CM, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); #View(merged_sigCM)

# # run biotype merging for each cell-type
# tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# merged_sig_biotypeCM = merge(x = merged_sigCM, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); #View(merged_sig_biotypeCM)


# #run merging with gene symbol
# tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# merged_sig_bio_symbCM = merge(x = merged_sig_biotypeCM, y = tx2symbol_mus, by.x = "tx", 
#                               by.y = "Transcript stable ID version", all.x = TRUE); View(merged_sig_bio_symbCM)


# write.table(merged_sig_bio_symbCM, file="sigprofiles_bio_sym_HC.txt", sep="\t", quote=F, col.names=NA)

# names(MyIso1$Tfit)
# View(MyIso1$Tfit$sig.profiles); dim(MyIso1$Tfit$sig.profiles)             #7621
# View(MyIso1$Tfit$group.coeffs); dim(MyIso1$Tfit$group.coeffs)         #7621

# View(MyIso1$Tfit$sig.profilesvariables)
# View(MyIso1$Tfit$sol); dim(MyIso1$Tfit$sol)

# # getDS creates lists of significant isoforms from Differentially Spliced Genes (DSG)
# Myget1 <- getDS(MyIso1)
# names(Myget1)
# names(Myget1$get2$sig.genes$Group5vsGroup0)
# names(Myget1$get2$summary)
# View(Myget1$get2$summary)


# see1 <- seeDS(Myget1, cluster.all=FALSE, k=6) #k.mclust=TRUE)
# Myget1$List0; Myget1$DET; Myget1$DSG
# # ============== plot plofiles of the clusters that differ in groups ==========================
# list1 = names(which(see1$cut==1))
# LIST1 <- data_ami[rownames(data_ami) %in% list1,-1]#dim(LIST1);View(LIST1)
# dev.off()
# PlotGroups(LIST1, edesign = dis_ami$edesign)
# PlotGroups(LIST1, edesign = dis_ami$edesign, show.fit = T, dis = dis_ami$dis, groups.vector = dis_ami$groups.vector)
# dev.off()
# PlotProfiles(LIST1, cond = rownames(edesign_ami), repvect = edesign_ami$replicate)
# ## ============================================================================================
# table <- tableDS(see1)

# PC<-PodiumChange(Myget1, only.sig.iso=TRUE, comparison="any"); #names(PC)
# PC$L; PC$edesign

# dev.off()
# IsoPlot(Myget1,"ENSMUSG00000073557.13",only.sig.iso=FALSE,cex.main=2,cex.legend=1)
# ## =================================== treating each transcript differnt ================================
# edesign_HC <- read.table(file.path("samples_HC_diff.txt"), header=TRUE)        # head(edesign_HC); 


# dis_ami <- make.design.matrix(edesign_HC, degree = 5)
# # p.vector performs a regression fit for each gene taking all variables present in the 
# # model given by a regression matrix and returns a list of FDR corrected significant genes.

# fit <- p.vector(data_ami[,-1], dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE);#names(fit)
# # fit <- p.vector(normcounts_HC, dis_ami, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts=TRUE);
# write.table(fit$SELEC, file="sigprofiles_bio_sym_HC_d14_d0.txt", sep="\t", quote=F, col.names=NA)

# tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)
# get<-get.siggenes(tstep, vars="groups")
# saveRDS(get, file = "MyIso_HC.rds")


# # get1 = readRDS(file = "MyIso_HC.rds") #get the RDS file from the server
get1 = readRDS(file = paste0("MyIso_gene_",opt$c,".rds")) #get the RDS file from the server
# # View(get1$summary)
# #======================================== save the results =================================================
# HC_d1_d0 = (cbind(get1$sig.genes$Group1vsGroup0$sig.profiles, get1$sig.genes$Group1vsGroup0$sig.pvalues$`p-value`))
# # HC_d1_d0$tx = rownames(HC_d1_d0)

# # run geneID merging for each cell-type
# # tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# # merged_diffHC_d1_d0 = merge(x = HC_d1_d0, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); 

# # run biotype merging for each cell-type
# # tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# # merged_diff_biotypeHC_d1_d0 = merge(x = merged_diffHC_d1_d0, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); 


# #run merging with gene symbol
# # tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# # merged_diff_bio_symbHC_d1_d0 = merge(x = merged_diff_biotypeHC_d1_d0, y = tx2symbol_mus, by.x = "tx", 
# #                               by.y = "Transcript stable ID version", all.x = TRUE); View(merged_diff_bio_symbHC_d1_d0)

# # write.table(merged_diff_bio_symbHC_d1_d0, file="diffprofiles_bio_sym_HC_d1_d0.txt", sep="\t", quote=F, col.names=NA)

# # for gene quantification 
# HC_d1_d0$gn = rownames(HC_d1_d0)
# #merge with gene symbol
# gene2biotype = read_delim(file.path("~/Downloads/mouse_gene_biotype.txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symbHC_d1_d0 = merge(x = HC_d1_d0, y = gene2biotype, by.x = "gn", 
#                                      by.y = "Gene stable ID version", all.x = TRUE); #View(merged_diff_bio_symbHC_d1_d0)
# write.table(merged_diff_bio_symbHC_d1_d0, file="diff_gene_profiles_bio_sym_HC_d1_d0.txt", sep="\t", quote=F, col.names=NA)
# #===============================================================================================================

# HC_d3_d0 = cbind(get1$sig.genes$Group2vsGroup0$sig.profiles, get1$sig.genes$Group2vsGroup0$sig.pvalues$`p-value`)
# # HC_d3_d0$tx = rownames(HC_d3_d0);
# # 
# # # run geneID merging for each cell-type
# # tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# # merged_diffHC_d3_d0 = merge(x = HC_d3_d0, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); 
# # 
# # # run biotype merging for each cell-type
# # tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# # merged_diff_biotypeHC_d3_d0 = merge(x = merged_diffHC_d3_d0, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); 
# # 
# # 
# # #run merging with gene symbol
# # tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# # merged_diff_bio_symbHC_d3_d0 = merge(x = merged_diff_biotypeHC_d3_d0, y = tx2symbol_mus, by.x = "tx", 
# #                                      by.y = "Transcript stable ID version", all.x = TRUE); View(merged_diff_bio_symbHC_d3_d0)
# # 
# # write.table(merged_diff_bio_symbHC_d3_d0, file="diffprofiles_bio_sym_HC_d3_d0.txt", sep="\t", quote=F, col.names=NA)

# # for gene quantification 
# HC_d3_d0$gn = rownames(HC_d3_d0)
# #merge with gene symbol
# gene2biotype = read_delim(file.path("~/Downloads/mouse_gene_biotype.txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symbHC_d3_d0 = merge(x = HC_d3_d0, y = gene2biotype, by.x = "gn", 
#                                      by.y = "Gene stable ID version", all.x = TRUE); #View(merged_diff_bio_symbHC_d3_d0)
# write.table(merged_diff_bio_symbHC_d3_d0, file="diff_gene_profiles_bio_sym_HC_d3_d0.txt", sep="\t", quote=F, col.names=NA)
# #===============================================================================================================
# HC_d7_d0 = cbind(get1$sig.genes$Group3vsGroup0$sig.profiles, get1$sig.genes$Group3vsGroup0$sig.pvalues$`p-value`)
# # HC_d7_d0$tx = rownames(HC_d7_d0);dim(HC_d7_d0)
# # 
# # # run geneID merging for each cell-type
# # tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# # merged_diffHC_d7_d0 = merge(x = HC_d7_d0, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); 
# # 
# # # run biotype merging for each cell-type
# # tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# # merged_diff_biotypeHC_d7_d0 = merge(x = merged_diffHC_d7_d0, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); 
# # 
# # 
# # #run merging with gene symbol
# # tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# # merged_diff_bio_symbHC_d7_d0 = merge(x = merged_diff_biotypeHC_d7_d0, y = tx2symbol_mus, by.x = "tx", 
# #                                      by.y = "Transcript stable ID version", all.x = TRUE); View(merged_diff_bio_symbHC_d7_d0)
# # 
# # write.table(merged_diff_bio_symbHC_d7_d0, file="diffprofiles_bio_sym_HC_d7_d0.txt", sep="\t", quote=F, col.names=NA)

# # for gene quantification 
# HC_d7_d0$gn = rownames(HC_d7_d0)
# #merge with gene symbol
# gene2biotype = read_delim(file.path("~/Downloads/mouse_gene_biotype.txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symbHC_d7_d0 = merge(x = HC_d7_d0, y = gene2biotype, by.x = "gn", 
#                                      by.y = "Gene stable ID version", all.x = TRUE); #View(merged_diff_bio_symbHC_d7_d0)
# write.table(merged_diff_bio_symbHC_d7_d0, file="diff_gene_profiles_bio_sym_HC_d7_d0.txt", sep="\t", quote=F, col.names=NA)
# #===============================================================================================================
# HC_d14_d0 = cbind(get1$sig.genes$Group4vsGroup0$sig.profiles, get1$sig.genes$Group4vsGroup0$sig.pvalues$`p-value`)
# # HC_d14_d0$tx = rownames(HC_d14_d0);dim(HC_d14_d0)
# # 
# # # run geneID merging for each cell-type
# # tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# # merged_diffHC_d14_d0 = merge(x = HC_d14_d0, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); 
# # 
# # # run biotype merging for each cell-type
# # tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# # merged_diff_biotypeHC_d14_d0 = merge(x = merged_diffHC_d14_d0, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); 
# # 
# # 
# # #run merging with gene symbol
# # tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# # merged_diff_bio_symbHC_d14_d0 = merge(x = merged_diff_biotypeHC_d14_d0, y = tx2symbol_mus, by.x = "tx", 
# #                                      by.y = "Transcript stable ID version", all.x = TRUE); View(merged_diff_bio_symbHC_d14_d0)
# # 
# # write.table(merged_diff_bio_symbHC_d14_d0, file="diffprofiles_bio_sym_HC_d14_d0.txt", sep="\t", quote=F, col.names=NA)

# # for gene quantification 
# HC_d14_d0$gn = rownames(HC_d14_d0)
# #merge with gene symbol
# gene2biotype = read_delim(file.path("~/Downloads/mouse_gene_biotype.txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symbHC_d14_d0 = merge(x = HC_d14_d0, y = gene2biotype, by.x = "gn", 
#                                      by.y = "Gene stable ID version", all.x = TRUE); #View(merged_diff_bio_symbHC_d14_d0)
# write.table(merged_diff_bio_symbHC_d14_d0, file="diff_gene_profiles_bio_sym_HC_d14_d0.txt", sep="\t", quote=F, col.names=NA)


# #===============================================================================================================
# HC_d28_d0 = cbind(get1$sig.genes$Group5vsGroup0$sig.profiles, get1$sig.genes$Group5vsGroup0$sig.pvalues$`p-value`)
# # HC_d28_d0$tx = rownames(HC_d28_d0);dim(HC_d28_d0)
# # 
# # # run geneID merging for each cell-type
# # tx2gene_mus = read_csv(file.path("tx2gene.csv"), show_col_types = FALSE)
# # merged_diffHC_d28_d0 = merge(x = HC_d28_d0, y = tx2gene_mus, by.x = "tx", by.y = "TXNAME", all.x = TRUE); 
# # 
# # # run biotype merging for each cell-type
# # tx2biotype_mus = read_csv(file.path("tx_biotype.csv"), show_col_types = FALSE)
# # merged_diff_biotypeHC_d28_d0 = merge(x = merged_diffHC_d28_d0, y = tx2biotype_mus, by.x = "tx", by.y = "tx", all.x = TRUE); 
# # 
# # 
# # #run merging with gene symbol
# # tx2symbol_mus = read_delim(file.path("~/Downloads/mart_export (1).txt"), show_col_types = TRUE, delim = "\t")
# # merged_diff_bio_symbHC_d28_d0 = merge(x = merged_diff_biotypeHC_d28_d0, y = tx2symbol_mus, by.x = "tx", 
# #                                       by.y = "Transcript stable ID version", all.x = TRUE); View(merged_diff_bio_symbHC_d28_d0)
# # 
# # write.table(merged_diff_bio_symbHC_d28_d0, file="diffprofiles_bio_sym_HC_d28_d0.txt", sep="\t", quote=F, col.names=NA)

# # for gene quantification 
# HC_d28_d0$gn = rownames(HC_d28_d0)
# #merge with gene symbol
# gene2biotype = read_delim(file.path("~/Downloads/mouse_gene_biotype.txt"), show_col_types = TRUE, delim = "\t")
# merged_diff_bio_symbHC_d28_d0 = merge(x = HC_d28_d0, y = gene2biotype, by.x = "gn", 
#                                      by.y = "Gene stable ID version", all.x = TRUE); #View(merged_diff_bio_symbHC_d28_d0)
# write.table(merged_diff_bio_symbHC_d28_d0, file="diff_gene_profiles_bio_sym_HC_d28_d0.txt", sep="\t", quote=F, col.names=NA)

# ======================================= install.packages("UpSetR") =============================================
library("UpSetR")
pdf(paste0("upset_",opt$c,".2.pdf"),
        width = 5, 
        height = 6,
        onefile=F)
listInput <- list(d1vsd0 = get1$summary$Group1vsGroup0, 
                  d3vsd0 = get1$summary$Group2vsGroup0, 
                  d7vsd0 = get1$summary$Group3vsGroup0,
                  d14vsd0 = get1$summary$Group4vsGroup0,
                  d28vsd0 = get1$summary$Group5vsGroup0)
upset(fromList(listInput), 
                order.by = "freq",
                text.scale = c(1.5, 1.5, 1.5, 0.7, 2, 0))
#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
dev.off()

## =================================== clusterProfiler on retained introns ======================================
#selected the retained intron genes , and made them as seperate files for each celltype
# setwd("/Users/ranjan/time_course_analysis/masigpro/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)

# ## load all gene symbols from each cluster profile
# profHC = read.delim(file="ri_HC", header = TRUE)
# profCM = read.delim(file="ri_CM", header = TRUE)
# profFB = read.delim(file="ri_FB", header = TRUE)
# profEC = read.delim(file="ri_EC", header = TRUE)

# ## convert to entrex gene symbol
# cp = list()
# cp$HC = bitr(profHC$ri_HC, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
# cp$CM = bitr(profCM$ri_CM, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
# cp$FB = bitr(profFB$ri_FB, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
# cp$EC = bitr(profEC$ri_EC, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID

# ck_go <- compareCluster(geneCluster = cp, fun = "enrichGO", OrgDb='org.Mm.eg.db' ) 
# dotplot(ck_go)


## ================= finding the common retained introns and NMD transcripts ==========================
#read the file (from Justin)



#read the transcript ids from cell type



# merge and count





