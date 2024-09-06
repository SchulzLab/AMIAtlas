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
write.table(normalized_counts, file="results/normalized_genecounts_HC.txt", sep="\t", quote=F, col.names=NA)


