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
# print("hello: Running miRNA differential expression analysis for EC samples")
paste("Current working directory:", getwd())

samples_EC <- read.table(file.path("EC/samples_EC_new.txt"), header=TRUE)
rownames(samples_EC) = c(paste("ctrl", rep("EC",4), 1:4, sep="_"),
  paste("d1", rep("EC",4), 1:4, sep="_"),
  paste("d3", rep("EC",4), 1:4, sep="_"),
  paste("d7", rep("EC",4), 1:4, sep="_"),
  paste("d14", rep("EC",4), 1:4, sep="_"),
  paste("d28", rep("EC",4), 1:4, sep="_"))


head(samples_EC,5)

data.mir.EC = read.table(file.path("EC/EC_mature_normalized_CPM.txt"), header=TRUE)

# head(data.mir.EC,5)
# design <- make.design.matrix(samples_EC, degree = 5)
# colnames(data.mir.EC) = rownames(samples_EC)
# data.mir.ECn = data.mir.EC[rowSums(data.mir.EC) > 10,]

# dim(data.mir.EC)
# dim(data.mir.ECn)

# fitEC <- p.vector(data = data.mir.ECn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

# fitEC$i # returns the number of significant genes
# paste("alpha", head(fitEC$alfa))# gives p-value at the Q false discovery control level
# head(fitEC$SELEC) # is a matrix with the significant genes and their expression values
# head(fitEC$p.vector) # vector containing the computed p-values
# head(fitEC$p.adjusted)

# names(fitEC)

# pheatmap(fitEC$SELEC,
#          scale="row",
#          clustering_distance_rows = "correlation",
#          cluster_cols=F,
#          show_rownames=F,
#          border_color=NA,
#          filename="Heatmap_significant_miR_EC_mod.png"
# )

# tEC <- T.fit(data = fitEC, step.method = "two.ways.forward")

# get<-get.siggenes(tEC, vars="groups")
# write.table(get$summary, "EC_summary_miRNAs.txt", na = "NA", quote=FALSE, sep=",")

# names(get$sig.genes)
# names(get$sig.genes$Group.1vsGroup.0)
# View(get$sig.genes$Group.1vsGroup.0$sig.pvalues)
# View(get$sig.genes$Group.1vsGroup.0$coefficients)
# names(get$summary$Group.1vsGroup.0)

# png("upset_EC_miR.png", width = 720, height = 720)
# listInput <- list(d1vsd0 = get$summary$Group.1vsGroup.0, 
#                   d3vsd0 = get$summary$Group.3vsGroup.0, 
#                   d7vsd0 = get$summary$Group.7vsGroup.0,
#                   d14vsd0 = get$summary$Group.14vsGroup.0,
#                   d28vsd0 = get$summary$Group.28vsGroup.0)
# upset(fromList(listInput), order.by = "freq")
# dev.off()
#==================== need to remove the outlier samples and run again ==================
print("hello: Running miRNA differential expression analysis removing outlier EC samples")
samples_EC <- read.table(file.path("EC/samples_EC_remout.txt"), header=TRUE)
rownames(samples_EC) = c(paste("ctrl", rep("EC",4), 1:4, sep="_"), 
                         paste("d1", rep("EC",4), 1:4, sep="_"),
                         paste("d3", rep("EC",4), 1:4, sep="_"),
                         paste("d7", rep("EC",4), 1:4, sep="_"),
                         paste("d14", rep("EC",4), 1:4, sep="_"),
                         paste("d28", rep("EC",2), 1:2, sep="_"),
                         paste("d28", rep("EC",1), 4, sep="_"))

head(samples_EC,5)

data.mir.EC = read.table(file.path("EC/EC_mature_normalized_CPM.txt"), header=TRUE)
drops <- c("d28_EC_3.mature")
data.mir.EC = data.mir.EC[ , !(colnames(data.mir.EC) %in% drops)]
head(data.mir.EC,3)

design <- make.design.matrix(samples_EC, degree = 5)
colnames(data.mir.EC) = rownames(samples_EC)
data.mir.ECn = data.mir.EC[rowSums(data.mir.EC) > 10,]

dim(data.mir.EC)
dim(data.mir.ECn)

fitEC <- p.vector(data = data.mir.ECn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

paste("# of sig. miRs", fitEC$i) # returns the number of significant genes
paste("alpha", head(fitEC$alfa))# gives p-value at the Q false discovery control level
head(fitEC$SELEC) # is a matrix with the significant genes and their expression values
head(fitEC$p.vector) # vector containing the computed p-values
head(fitEC$p.adjusted)

pheatmap(fitEC$SELEC,
         scale="row",
         clustering_distance_rows = "correlation",
         cluster_cols=F,
         show_rownames=F,
         border_color=NA,
         #filename="Heatmap_significant_miR_EC_remout.png"
)

tEC <- T.fit(data = fitEC, step.method = "two.ways.forward")

get<-get.siggenes(tEC, vars="groups")
write.table(get$summary, "EC_summary_miRNAs_remout.txt", na = "NA", quote=FALSE, sep=",")

# get$sig.genes$Group.1vsGroup.0$g
# get$sig.genes$Group.3vsGroup.0$g
# get$sig.genes$Group.7vsGroup.0$g
# get$sig.genes$Group.14vsGroup.0$g
# get$sig.genes$Group.28vsGroup.0$g

# names(get$sig.genes)
# names(get$sig.genes$Group.1vsGroup.0)
# View(get$sig.genes$Group.1vsGroup.0$sig.pvalues)
# View(get$sig.genes$Group.1vsGroup.0$coefficients)
# names(get$summary$Group.1vsGroup.0)

pdf("upset_EC_miR_remout.1.pdf", width = 4, height = 6, onefile=F)
listInput <- list(d1vsd0 = get$summary$Group.1vsGroup.0, 
                  d3vsd0 = get$summary$Group.3vsGroup.0, 
                  d7vsd0 = get$summary$Group.7vsGroup.0,
                  d14vsd0 = get$summary$Group.14vsGroup.0,
                  d28vsd0 = get$summary$Group.28vsGroup.0)
upset(fromList(listInput), 
      order.by = "freq", 
      sets.x.label = "DE miRNAs per timepoint",
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 0))
dev.off()
#========================================================================================