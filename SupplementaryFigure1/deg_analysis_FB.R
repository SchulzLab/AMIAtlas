library(tximeta)
library(tximport)
suppressWarnings(suppressMessages(library("DESeq2")))
library(ggplot2)

dir = getwd()

# samples_FB <- read.table(file.path("samples_FB.txt"), header=TRUE)
samples_FB <- read.table(file.path("samples_FB.new.txt"), header=TRUE)        #<-------------- removed outlier samples

filesFB <- file.path(dir, samples_FB$pop, paste(samples_FB$run, "_transcripts_quant", sep=""), "quant.sf")

names(filesFB) = c(paste("d0", rep("FB",4), 1:4, sep="_"),
               paste("d1", rep("FB",4), 1:4, sep="_"),
               # paste("d3", rep("FB",4), 1:4,sep="_"), 
               paste("d3", rep("FB",2), c(1,4),sep="_"),              #<-------------- removed outlier samples            
               paste("d7", rep("FB",4), 1:4, sep="_"),
               paste("d14", rep("FB",4), 1:4, sep="_"),
               # paste("d14", rep("FB",4), paste(1:4, ".", sep="", 2), sep="_"),
               paste("d28", rep("FB",4), 1:4, sep="_"))
               # paste("d28", rep("FB",4), paste(1:4, ".", sep="", 2), sep="_"))

txi <- tximport(filesFB, type="salmon", txOut=TRUE) 


## =================== on gene level =========================
library(readr)
tx2gene_mus = read_csv(file.path("tx2gene.csv"))

# head(tx2gene_mus) 
txigene <- tximport(filesFB, type="salmon", tx2gene=tx2gene_mus)

samples_FB$timepoint <- as.factor(samples_FB$timepoint)
dds <- DESeqDataSetFromTximport(txi = txigene,           #txi, 
                              colData = samples_FB,
                              design = ~ rep + timepoint)

# dds
paste("before count filtering:", nrow(dds))
# # at least 3 samples with a count of 10 or higher
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

# ggsave("SampleDistances_FB_genelevel_woutoutliers.png", plot = p)
## ================================ PCA plots =====================================
# head(counts(dds))
cat("plotting PCA remofing outliers ...")
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("timepoint", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, shape=timepoint, color=rep)) +
  geom_point(size=3) +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("FB: PCA with VST data (on Gene Level)")

ggsave("PCAplot_FB_genelevel_woutoutliers.pdf", plot = pca)

## ===============normalization : generate normalized counts, then TPM =============
cat("perform DESEQ2 normalization ....")
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)                       # summary(normalized_counts)
normalized_counts_df = as.data.frame(normalized_counts)
normalized_counts_df$gname = rownames(normalized_counts_df)            # head(normalized_counts_df)

# colnames(normalized_counts_df)
write.table(normalized_counts, file="data/normalized_genecounts_FB.txt", sep="\t", quote=F, col.names=NA)
