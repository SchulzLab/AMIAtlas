library(tximeta)
library(tximport)
suppressWarnings(suppressMessages(library("DESeq2")))
library(ggplot2)

dir = getwd()

samples_CM <- read.table(file.path("samples_CM.txt"), header=TRUE)

filesCM <- file.path(dir, samples_CM$pop, paste(samples_CM$run, "_transcripts_quant", sep=""), "quant.sf")

names(filesCM) = c(paste("d0", rep("CM",4), 1:4, sep="_"),
               paste("d1", rep("CM",4), 1:4, sep="_"),
               paste("d3", rep("CM",4), 1:4, sep="_"),
               # paste("d3", rep("CM",3), paste(4, ".", sep="", 1:3), sep="_"),
               paste("d7", rep("CM",4), 1:4, sep="_"),
               paste("d14", rep("CM",4), 1:4, sep="_"),
               # paste("d14", rep("CM",4), paste(1:4, ".", sep="", 2), sep="_"),
               paste("d28", rep("CM",4), 1:4, sep="_"))
               # paste("d28", rep("CM",4), paste(1:4, ".", sep="", 2), sep="_"))

# txi <- tximport(filesCM, type="salmon", txOut=TRUE) 

## =================== on gene/transcript level =========================
library(readr)
tx2gene_mus = read_csv(file.path("tx2gene.csv"))
txigene <- tximport(filesCM, type="salmon", tx2gene=tx2gene_mus)


samples_CM$timepoint <- as.factor(samples_CM$timepoint)
dds <- DESeqDataSetFromTximport(txi = txigene,                      #txi,
                              colData = samples_CM,
                              design = ~ rep + timepoint)

# filtering
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
# rownames(samplePoisDistMatrix) <- paste( dds$timepoint, dds$rep, sep=" - " )
# colnames(samplePoisDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# p = pheatmap(samplePoisDistMatrix,
#          clustering_distance_rows = poisd$dd,
#          clustering_distance_cols = poisd$dd,
#          col = colors)

# ggsave("SampleDistances_CM_genelevel.png", plot = p)

## ========================== PCA plots  ===============================
head(counts(dds))

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("timepoint", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, shape=timepoint, color=rep)) +
  geom_point(size=3) +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("CM: PCA with VST data")

ggsave("PCAplot_CM.pdf", plot = pca)


## ===============normalization : generate normalized counts =============
print("perform DESEQ2 normalization ....")
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)                       # summary(normalized_counts)
normalized_counts_df = as.data.frame(normalized_counts)
normalized_counts_df$txname = rownames(normalized_counts_df)            # head(normalized_counts_df)

# colnames(normalized_counts_df)
write.table(normalized_counts, file="data/normalized_genecounts_CM.txt", sep="\t", quote=F, col.names=NA)

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

# c_addgene = c();medgene = c(); sd_gene = c();

# c_addpcg = c(); c_addpse = c(); c_addri = c(); c_addlnc = c(); c_addnmd = c(); c_addpt = c(); c_addsn = c(); c_addsno = c(); c_addr = c(); c_addmi = c()
# medpcg = c(); medpse = c(); medri = c(); medlnc = c(); mednmd = c(); medpt = c(); medsn = c(); medsno = c(); medr = c(); medmi = c()
# sd_pcg = c(); sd_pse = c(); sd_ri = c(); sd_lnc = c(); sd_nmd = c(); sd_pt = c(); sd_sn = c(); sd_sno = c(); sd_r = c(); sd_mi = c()

# for (i in 3:length(colnames(merged))) {

#     count = count+1

#     if(count %% 4 == 0){  # calculate median after every 4th sample

#       #add the last element to the list before median calculation
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
      
#       print(paste0("compute median here: pcg", count ,"median: ", median(c_addpcg), "Standard Dev:", sd(c_addpcg)))    # print(c_addpcg)
#       print(paste0("compute median here: pse", count ,"median: ", median(c_addpse), "Standard Dev:", sd(c_addpse)))
#       print(paste0("compute median here: ri", count ,"median: ", median(c_addri), "Standard Dev:", sd(c_addri)))
#       print(paste0("compute median here: lncRNA", count ,"median: ", median(c_addlnc), "Standard Dev:", sd(c_addlnc)))
#       print(paste0("compute median here: nmd", count ,"median: ", median(c_addnmd), "Standard Dev:", sd(c_addnmd)))
#       print(paste0("compute median here: pt", count ,"median: ", median(c_addpt), "Standard Dev:", sd(c_addpt)))
#       print(paste0("compute median here: snRNA", count ,"median: ", median(c_addsn), "Standard Dev:", sd(c_addsn)))
#       print(paste0("compute median here: snoRNA", count ,"median: ", median(c_addsno), "Standard Dev:", sd(c_addsno)))
#       print(paste0("compute median here: rRNA", count ,"median: ", median(c_addr), "Standard Dev:", sd(c_addr)))
#       print(paste0("compute median here: miRNA", count ,"median: ", median(c_addmi), "Standard Dev:", sd(c_addmi)))


#       medpcg = c(medpcg, median(c_addpcg)); sd_pcg = c(sd_pcg, sd(c_addpcg)); 
#       medpse = c(medpse, median(c_addpse)); sd_pse = c(sd_pse, sd(c_addpse));
#       medri = c(medri, median(c_addri)); sd_ri = c(sd_ri, sd(c_addri));
#       medlnc = c(medlnc, median(c_addlnc)); sd_lnc = c(sd_lnc, sd(c_addlnc));
#       mednmd = c(mednmd, median(c_addnmd)); sd_nmd = c(sd_nmd, sd(c_addnmd));
#       medpt = c(medpt, median(c_addpt)); sd_pt = c(sd_pt, sd(c_addpt));
#       medsn = c(medsn, median(c_addsn)); sd_sn = c(sd_sn, sd(c_addsn));
#       medsno = c(medsno, median(c_addsno)); sd_sno = c(sd_sno, sd(c_addsno));
#       medr = c(medr, median(c_addr)); sd_r = c(sd_r, sd(c_addr));
#       medmi = c(medmi, median(c_addmi)); sd_mi = c(sd_mi, sd(c_addmi));

#       medgene = c(medgene, median(c_addgene)); sd_gene = c(sd_gene, sd(c_addgene));

#       c_addpcg = c(); c_addpse = c(); c_addri = c(); c_addlnc = c(); c_addnmd = c(); c_addpt = c(); 
#       c_addsn = c(); c_addsno = c(); c_addr = c(); c_addmi = c(); c_addgene = c()

#     }else{      # append to median calculating list

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
#     }


# }


# merged_calc = data.frame(medpcg, sd_pcg, medpse, sd_pse, medri, sd_ri, medlnc, sd_lnc, mednmd, sd_nmd, medpt, sd_pt, medsn,sd_sn, medsno, sd_sno, medr, sd_r, medmi, sd_mi, medgene, sd_gene) 

# names(merged_calc) <- c('medpcg', 'sd_pcg', 'medpse', 'sd_pse', 'medri', 'sd_ri', 'medlnc', 'sd_lnc', 'mednmd', 
#   'sd_nmd', 'medpt', 'sd_pt', 'medsn','sd_sn', 'medsno', 'sd_sno', 'medr', 'sd_r', 'medmi', 'sd_mi', 'medgene', 'sd_gene')

# merged_calc
# write.table(merged, file="data/normalized_counts_CM_biotype.txt", sep="\t", quote=F, col.names=NA)


## ================================ Violin plots ======================================

# suppressWarnings(suppressMessages(library(tidyr)))
# suppressWarnings(suppressMessages(library(ggplot2)))
# suppressWarnings(suppressMessages(library(dplyr)))
# suppressWarnings(suppressMessages(library(ggpubr)))

# data_wide <- as.data.frame(as.matrix(assay(vsd))); rownames(data_wide) <- NULL

# violin_d0 <- data_wide[,1:4] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D0") +
#     geom_violin()

# violin_d1 <- data_wide[,5:8] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D1") +
#     geom_violin()
 
# ggsave("Violinplot_CM_d1.png", plot = violin_d1)

# violin_d3 <- data_wide[,9:12] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D3") +
#     geom_violin()
 
# ggsave("Violinplot_CM_d3.png", plot = violin_d3)

# violin_d7 <- data_wide[,13:16] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D7") +
#     geom_violin()
 
# ggsave("Violinplot_CM_d7.png", plot = violin_d7)

# violin_d14 <- data_wide[,17:20] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D14") +
#     geom_violin()
 
# ggsave("Violinplot_CM_d14.png", plot = violin_d14)

# violin_d28 <- data_wide[,21:24] %>% 
#   gather(key="Replicates", value="Val") %>%
#   ggplot( aes(x=Replicates, y=Val, fill=Replicates)) +
#   ggtitle("D28") +
#     geom_violin()
 
# ggsave("Violinplot_CM_d28.png", plot = violin_d28)

# save = ggarrange(violin_d0,violin_d1,violin_d3, violin_d7, violin_d14, violin_d28,
#           labels = c("i", "ii", "iii", "iv", "v", "vi"),
#           ncol = 6, nrow = 1)
# ggsave("Violinplot_CM.png", plot = save)