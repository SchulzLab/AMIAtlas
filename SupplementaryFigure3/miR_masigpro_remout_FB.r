#!/usr/bin/env Rscript
#========================== package installation scripts ==============================
if (!require("maSigPro")){
	BiocManager::install("maSigPro")

	library(maSigPro)
}
if (!require("UpSetR")){
	install.packages("UpSetR")
	library("UpSetR")
}
library(pheatmap)

#=================================================================================
setwd("/projects/amitimeseries/work/smallRNA_AMI/miRBase_counts/")
# print("hello: Running miRNA differential expression analysis for FB samples")
paste("Current working directory:", getwd())

samples_FB <- read.table(file.path("HC/samples_HC_new.txt"), header=TRUE)
rownames(samples_FB) = c(paste("ctrl", rep("FB",4), 1:4, sep="_"),
  paste("d1", rep("FB",4), 1:4, sep="_"),
  paste("d3", rep("FB",4), 1:4, sep="_"),
  paste("d7", rep("FB",4), 1:4, sep="_"),
  paste("d14", rep("FB",4), 1:4, sep="_"),
  paste("d28", rep("FB",4), 1:4, sep="_"))


# head(samples_FB,5)

# data.mir.FB = read.table(file.path("FB/FB_mature_normalized_CPM.txt"), header=TRUE)

# head(data.mir.FB,5)
# design <- make.design.matrix(samples_FB, degree = 5)
# colnames(data.mir.FB) = rownames(samples_FB)
# data.mir.FBn = data.mir.FB[rowSums(data.mir.FB) > 10,]

# dim(data.mir.FB)
# dim(data.mir.FBn)

# fitFB <- p.vector(data = data.mir.FBn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

# fitFB$i # returns the number of significant genes
# paste("alpha", head(fitFB$alfa))# gives p-value at the Q false discovery control level
# head(fitFB$SELEC) # is a matrix with the significant genes and their expression values
# head(fitFB$p.vector) # vector containing the computed p-values
# head(fitFB$p.adjusted)

# names(fitFB)

# pheatmap(fitFB$SELEC,
#          scale="row",
#          clustering_distance_rows = "correlation",
#          cluster_cols=F,
#          show_rownames=F,
#          border_color=NA,
#          filename="Heatmap_significant_miR_FB_mod.png"
# )

# tFB <- T.fit(data = fitFB, step.method = "two.ways.forward")

# get<-get.siggenes(tFB, vars="groups")
# write.table(get$summary, "FB_summary_miRNAs.txt", na = "NA", quote=FALSE, sep=",")

# names(get$sig.genes)
# names(get$sig.genes$Group.1vsGroup.0)
# View(get$sig.genes$Group.1vsGroup.0$sig.pvalues)
# View(get$sig.genes$Group.1vsGroup.0$coefficients)
# names(get$summary$Group.1vsGroup.0)

# png("upset_FB_miR.png", width = 600, height = 600, res=100)
# listInput <- list(d1vsd0 = get$summary$Group.1vsGroup.0, 
#                   d3vsd0 = get$summary$Group.3vsGroup.0, 
#                   d7vsd0 = get$summary$Group.7vsGroup.0,
#                   d14vsd0 = get$summary$Group.14vsGroup.0,
#                   d28vsd0 = get$summary$Group.28vsGroup.0)
# upset(fromList(listInput), order.by = "freq")
# dev.off()
#==================== need to remove the outlier samples and run again d0, rep 2 ==================
print("hello: Running miRNA differential expression analysis removing outlier FB samples")
samples_FB <- read.table(file.path("FB/samples_FB_remout.txt"), header=TRUE)
rownames(samples_FB) = c(paste("ctrl", rep("FB",1), 1, sep="_"), 
                         paste("ctrl", rep("FB",2), 3:4, sep="_"),
                         paste("d1", rep("FB",4), 1:4, sep="_"),
                         paste("d3", rep("FB",4), 1:4, sep="_"),
                         paste("d7", rep("FB",4), 1:4, sep="_"),
                         paste("d14", rep("FB",4), 1:4, sep="_"),
                         paste("d28", rep("FB",4), 1:4, sep="_"))

head(samples_FB,5)

data.mir.FB = read.table(file.path("FB/FB_mature_normalized_CPM.txt"), header=TRUE)
drops <- c("ctrl_FB_2.mature")
data.mir.FB = data.mir.FB[ , !(colnames(data.mir.FB) %in% drops)]
head(data.mir.FB,3)

design <- make.design.matrix(samples_FB, degree = 5)
colnames(data.mir.FB) = rownames(samples_FB)
data.mir.FBn = data.mir.FB[rowSums(data.mir.FB) > 10,]

dim(data.mir.FB)
dim(data.mir.FBn)

fitFB <- p.vector(data = data.mir.FBn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

paste("# of sig. miRs", fitFB$i) # returns the number of significant genes
paste("alpha", head(fitFB$alfa))# gives p-value at the Q false discovery control level
head(fitFB$SELEC) # is a matrix with the significant genes and their expression values
head(fitFB$p.vector) # vector containing the computed p-values
head(fitFB$p.adjusted)

pheatmap(fitFB$SELEC,
         scale="row",
         clustering_distance_rows = "correlation",
         cluster_cols=F,
         show_rownames=F,
         border_color=NA,
         filename="Heatmap_significant_miR_FB_remout.png"
)

tFB <- T.fit(data = fitFB, step.method = "two.ways.forward")

get<-get.siggenes(tFB, vars="groups")
write.table(get$summary, "FB_summary_miRNAs_remout.txt", na = "NA", quote=FALSE, sep=",")

# names(get$sig.genes)
# names(get$sig.genes$Group.1vsGroup.0)
# View(get$sig.genes$Group.1vsGroup.0$sig.pvalues)
# View(get$sig.genes$Group.1vsGroup.0$coefficients)
# names(get$summary$Group.1vsGroup.0)

pdf("upset_FB_miR_remout.1.pdf", 
    width = 4, height = 6, onefile=F)
listInput <- list(d1vsd0 = get$summary$Group.1vsGroup.0, 
                  d3vsd0 = get$summary$Group.3vsGroup.0, 
                  d7vsd0 = get$summary$Group.7vsGroup.0,
                  d14vsd0 = get$summary$Group.14vsGroup.0,
                  d28vsd0 = get$summary$Group.28vsGroup.0)
upset(fromList(listInput), order.by = "freq", 
      sets.x.label = "DE miRNAs per timepoint",
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 0))
dev.off()
#========================================================================================