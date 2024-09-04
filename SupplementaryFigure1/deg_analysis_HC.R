library(tximeta)
library(tximport)
suppressWarnings(suppressMessages(library("DESeq2")))
library(ggplot2)

dir = getwd()

# samples_HC <- read.table(file.path("samples_HC.txt"), header=TRUE)
samples_HC <- read.table(file.path("samples_HC.new.txt"), header=TRUE)

filesHC <- file.path(dir, samples_HC$pop, paste(samples_HC$run, "_transcripts_quant", sep=""), "quant.sf")

names(filesHC) = c(paste("d0", rep("HC",4), 1:4, sep="_"),
               paste("d1", rep("HC",4), 1:4, sep="_"),
               # paste("d3", rep("HC",4), 1:4, sep="_"),
               paste("d3", rep("HC",3), 1:3, sep="_"),
               paste("d7", rep("HC",4), 1:4, sep="_"),
               paste("d14", rep("HC",4), 1:4, sep="_"),
               paste("d28", rep("HC",3), 1:3, sep="_")
               # paste("d28", rep("HC",4), 1:4, sep="_")
               )
               
txi <- tximport(filesHC, type="salmon", txOut=TRUE) 

## ========================== on gene level ========================
library(readr)
tx2gene_mus = read_csv(file.path("tx2gene.csv"))

txigene <- tximport(filesHC, type="salmon", tx2gene=tx2gene_mus)

samples_HC$timepoint <- as.factor(samples_HC$timepoint)
dds <- DESeqDataSetFromTximport(txi = txigene,             #txi,
                              colData = samples_HC,
                              design = ~ rep + timepoint)

# dds
paste("before count filtering:", nrow(dds))
# at least 3 samples with a count of 10 or higher
# keep <- rowSums(counts(dds) >= 10) >= 3
# dds <- dds[keep,]
# paste("after count filtering:", nrow(dds))

## ========================== sample distances ========================
# library("PoiClaClu")
# library("pheatmap")
# library("RColorBrewer")
# poisd <- PoissonDistance(t(counts(dds)))

# samplePoisDistMatrix <- as.matrix( poisd$dd )
# rownames(samplePoisDistMatrix) <- paste( paste("D",dds$timepoint), paste("rep",dds$rep), sep=" - " )
# colnames(samplePoisDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# p = pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = colors)


# ggsave("SampleDistances_HC_genelevel_woutoutliers.png", plot = p)
## ========================== PCA plots (on transcript level) ===============================	

head(counts(dds))
nsub=nrow(dds)
vsd <- vst(dds, blind=FALSE)      #varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("timepoint", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, shape=timepoint, color=rep)) +
  geom_point(size=3) +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("HC: PCA with VST data (on gene level)")

ggsave("PCAplot_HC_genelevel_woutoutliers.pdf", plot = pca)

## ===============normalization : generate normalized counts =============
cat("perform DESEQ2 normalization ....")
dds <- estimateSizeFactors(dds)

# sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)                       # summary(normalized_counts)

normalized_counts_df = as.data.frame(normalized_counts)

normalized_counts_df$txname = rownames(normalized_counts_df)            # head(normalized_counts_df)

colnames(normalized_counts_df)
write.table(normalized_counts, file="data/normalized_genecounts_HC.txt", sep="\t", quote=F, col.names=NA)
# write.table(normalized_counts, file="data/normalized_counts_HC.txt", sep="\t", quote=F, col.names=NA)


# counts.mat = as.matrix(normalized_counts)

# gene.length = as.matrix(txi$length)

# x <- counts.mat / gene.length

# tpm.mat <- t( t(x) * 1e6 / colSums(x) )

# paste("--------------------------head of TPM matrix-------------------------")
# head(tpm.mat)

# paste("--------------------------summary of TPM matrix-------------------------")
# summary(tpm.mat)

# tpm.mat_df = as.data.frame(tpm.mat)
# tpm.mat_df$txname = rownames(normalized_counts_df);#head(tpm.mat_df)

## ================= map to transcript annotation file and count ===================

# tx2biotype_mus = read_csv(file.path("tx_biotype.csv"))              #mapping file: tx and annotation

# merged = merge(x=tx2biotype_mus, y=tpm.mat_df, by.x = "tx", by.y = "txname")
# merged_calc = data.frame()
# # # write.table(merged, file="data/normalized_counts_CM_biotype.txt", sep="\t", quote=F, col.names=NA)

# library(stringr)

# thrsh = 1; count=0;
# addpcg = c(); addpse = c(); addri = c(); addlnc = c(); addnmd = c(); addpt = c(); addsn = c(); addsno = c(); addr = c(); addmi = c()

# c_addgene = c();medgene = c(); sd_gene = c();
# c_addpcg = c(); c_addpse = c(); c_addri = c(); c_addlnc = c(); c_addnmd = c(); c_addpt = c(); c_addsn = c(); c_addsno = c(); c_addr = c(); c_addmi = c()
# medpcg = c(); medpse = c(); medri = c(); medlnc = c(); mednmd = c(); medpt = c(); medsn = c(); medsno = c(); medr = c(); medmi = c()
# sd_pcg = c(); sd_pse = c(); sd_ri = c(); sd_lnc = c(); sd_nmd = c(); sd_pt = c(); sd_sn = c(); sd_sno = c(); sd_r = c(); sd_mi = c()

# for (i in 3:length(colnames(merged))) {

#     count = count+1

#     # if(count %% 4 == 0){  # calculate 

#     #   #add the last element to the list before median calculation
#     #   c_addpcg = c(c_addpcg, (sum(merged$biotype == "protein_coding" & merged[[i]] > thrsh)/58823) * 100)
#     #   c_addpse = c(c_addpse, (sum(str_detect(merged$biotype, "pseudogene") & merged[[i]] > thrsh)/13479) * 100)
#     #   c_addri = c(c_addri, (sum(merged$biotype == "retained_intron" & merged[[i]] > thrsh)/21920) * 100)
#     #   c_addlnc = c(c_addlnc, (sum(merged$biotype == "lncRNA" & merged[[i]] > thrsh)/14716) * 100)
#     #   c_addnmd = c(c_addnmd, (sum(merged$biotype == "nonsense_mediated_decay" & merged[[i]] > thrsh)/7197) * 100)
#     #   c_addpt = c(c_addpt, (sum(merged$biotype == "processed_transcript" & merged[[i]] > thrsh)/15032) * 100)
#     #   c_addsn = c(c_addsn, (sum(merged$biotype == "snRNA" & merged[[i]] > thrsh)/1305) * 100)
#     #   c_addsno = c(c_addsno, (sum(merged$biotype == "snoRNA" & merged[[i]] > thrsh)/1364) * 100)
#     #   c_addr = c(c_addr, (sum(merged$biotype == "rRNA" & merged[[i]] > thrsh)/306) * 100)
#     #   c_addmi = c(c_addmi, (sum(merged$biotype == "miRNA" & merged[[i]] > thrsh)/2022) * 100)

#     #   c_addgene = c(c_addgene, (sum(str_detect(merged$biotype, "_gene") & merged[[i]] > thrsh)/626) * 100)
      
#     #   print(paste0("compute median here: pcg", count ,"median: ", median(c_addpcg), "Standard Dev:", sd(c_addpcg)))    # print(c_addpcg)
#     #   print(paste0("compute median here: pse", count ,"median: ", median(c_addpse), "Standard Dev:", sd(c_addpse)))
#     #   print(paste0("compute median here: ri", count ,"median: ", median(c_addri), "Standard Dev:", sd(c_addri)))
#     #   print(paste0("compute median here: lncRNA", count ,"median: ", median(c_addlnc), "Standard Dev:", sd(c_addlnc)))
#     #   print(paste0("compute median here: nmd", count ,"median: ", median(c_addnmd), "Standard Dev:", sd(c_addnmd)))
#     #   print(paste0("compute median here: pt", count ,"median: ", median(c_addpt), "Standard Dev:", sd(c_addpt)))
#     #   print(paste0("compute median here: snRNA", count ,"median: ", median(c_addsn), "Standard Dev:", sd(c_addsn)))
#     #   print(paste0("compute median here: snoRNA", count ,"median: ", median(c_addsno), "Standard Dev:", sd(c_addsno)))
#     #   print(paste0("compute median here: rRNA", count ,"median: ", median(c_addr), "Standard Dev:", sd(c_addr)))
#     #   print(paste0("compute median here: miRNA", count ,"median: ", median(c_addmi), "Standard Dev:", sd(c_addmi)))


#     #   medpcg = c(medpcg, median(c_addpcg)); sd_pcg = c(sd_pcg, sd(c_addpcg)); 
#     #   medpse = c(medpse, median(c_addpse)); sd_pse = c(sd_pse, sd(c_addpse));
#     #   medri = c(medri, median(c_addri)); sd_ri = c(sd_ri, sd(c_addri));
#     #   medlnc = c(medlnc, median(c_addlnc)); sd_lnc = c(sd_lnc, sd(c_addlnc));
#     #   mednmd = c(mednmd, median(c_addnmd)); sd_nmd = c(sd_nmd, sd(c_addnmd));
#     #   medpt = c(medpt, median(c_addpt)); sd_pt = c(sd_pt, sd(c_addpt));
#     #   medsn = c(medsn, median(c_addsn)); sd_sn = c(sd_sn, sd(c_addsn));
#     #   medsno = c(medsno, median(c_addsno)); sd_sno = c(sd_sno, sd(c_addsno));
#     #   medr = c(medr, median(c_addr)); sd_r = c(sd_r, sd(c_addr));
#     #   medmi = c(medmi, median(c_addmi)); sd_mi = c(sd_mi, sd(c_addmi));

#     #   medgene = c(medgene, median(c_addgene)); sd_gene = c(sd_gene, sd(c_addgene));

#     #   c_addpcg = c(); c_addpse = c(); c_addri = c(); c_addlnc = c(); c_addnmd = c(); c_addpt = c(); 
#     #   c_addsn = c(); c_addsno = c(); c_addr = c(); c_addmi = c(); c_addgene = c()

#     # }else{      # append to median calculating list

#       c_addpcg = c(c_addpcg, (sum(merged$biotype == "protein_coding" & merged[[i]] > thrsh)/58823) * 100)
#       c_addpse = c(c_addpse, (sum(str_detect(merged$biotype, "pseudogene") & merged[[i]] > thrsh)/13479) * 100)
#       c_addri = c(c_addri, (sum(merged$biotype == "retained_intron" & merged[[i]] > thrsh)/21920) * 100)
#       c_addlnc = c(c_addlnc, (sum(merged$biotype == "lncRNA" & merged[[i]] > thrsh)/14716) * 100)
#       c_addnmd = c(c_addnmd, (sum(merged$biotype == "nonsense_mediated_decay" & merged[[i]] > thrsh)/7197) * 100)
#       c_addpt = c(c_addpt, (sum(merged$biotype == "processed_transcript" & merged[[i]] > thrsh)/15032) * 100)
#       c_addsn = c(c_addsn, (sum(merged$biotype == "snRNA" & merged[[i]] > thrsh)/1305) * 100)
#       c_addsno = c(c_addsno, (sum(merged$biotype == "snoRNA" & merged[[i]] > thrsh)/1364) * 100)
#       c_addr = c(c_addr, (sum(merged$biotype == "rRNA" & merged[[i]] > thrsh)/306) * 100)
#       c_addmi = c(c_addmi, (sum(merged$biotype == "miRNA" & merged[[i]] > thrsh)/2022) * 100)

#       c_addgene = c(c_addgene, (sum(str_detect(merged$biotype, "_gene") & merged[[i]] > thrsh)/626) * 100)
#     # }
    

# }


# merged_calc = data.frame(c_addpcg,c_addpse,c_addri ,c_addlnc ,c_addnmd ,c_addpt ,c_addsn ,c_addsno ,c_addr ,c_addmi ,c_addgene)

# # merged_calc = data.frame(medpcg, sd_pcg, medpse, sd_pse, medri, sd_ri, medlnc, sd_lnc, mednmd, sd_nmd, medpt, sd_pt, 
# #   medsn,sd_sn, medsno, sd_sno, medr, sd_r, medmi, sd_mi, medgene, sd_gene) 

# names(merged_calc) <- c('c_addpcg','c_addpse','c_addri' ,'c_addlnc' ,'c_addnmd' ,'c_addpt' ,'c_addsn' ,'c_addsno' ,'c_addr' ,'c_addmi', 'c_addgene')

# # names(merged_calc) <- c('medpcg', 'sd_pcg', 'medpse', 'sd_pse', 'medri', 'sd_ri', 'medlnc', 'sd_lnc', 'mednmd', 
# #   'sd_nmd', 'medpt', 'sd_pt', 'medsn','sd_sn', 'medsno', 'sd_sno', 'medr', 'sd_r', 'medmi', 'sd_mi', 'medgene', 'sd_gene')

# merged_calc

# colnames(merged)

## ========================== Violin plots ===============================

# suppressWarnings(suppressMessages(library(tidyr)))
# suppressWarnings(suppressMessages(library(ggplot2)))
# suppressWarnings(suppressMessages(library(dplyr)))


# data_wide <- as.data.frame(as.matrix(assay(vsd))); rownames(data_wide) <- NULL

# violin <- data_wide[,1:4] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D0") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d0_immuno.png", plot = violin)

# violin <- data_wide[,5:8] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D1") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d1_immuno.png", plot = violin)

# violin <- data_wide[,9:12] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D3") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d3_immuno.png", plot = violin)

# violin <- data_wide[,13:16] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D7") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d7_immuno.png", plot = violin)

# violin <- data_wide[,17:20] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D14") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d14_immuno.png", plot = violin)

# violin <- data_wide[,21:24] %>% 
#   gather(key="MesureType", value="Val") %>%
#   ggplot( aes(x=MesureType, y=Val, fill=MesureType)) + ggtitle("D28") +
#     geom_violin()
 
# ggsave("Violinplot_HC_d28_immuno.png", plot = violin)