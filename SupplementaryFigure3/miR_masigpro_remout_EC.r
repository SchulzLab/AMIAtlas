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

#==================== remove the outlier samples and run again ==================
print("hello: Running miRNA differential expression analysis removing outlier EC samples")
samples_EC <- read.table(file.path("./configurations/samples_EC_remout.txt"), header=TRUE)
rownames(samples_EC) = c(paste("ctrl", rep("EC",4), 1:4, sep="_"), 
                         paste("d1", rep("EC",4), 1:4, sep="_"),
                         paste("d3", rep("EC",4), 1:4, sep="_"),
                         paste("d7", rep("EC",4), 1:4, sep="_"),
                         paste("d14", rep("EC",4), 1:4, sep="_"),
                         paste("d28", rep("EC",2), 1:2, sep="_"),
                         paste("d28", rep("EC",1), 4, sep="_"))

data.mir.EC = read.table(file.path("./expression/EC/EC_mature_normalized_CPM.txt"), header=TRUE)
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

pheatmap(fitEC$SELEC,
         scale="row",
         clustering_distance_rows = "correlation",
         cluster_cols=F,
         show_rownames=F,
         border_color=NA,
         filename="Heatmap_significant_miR_EC_remout.png"
)

tEC <- T.fit(data = fitEC, step.method = "two.ways.forward")

get<-get.siggenes(tEC, vars="groups")
write.table(get$summary, "EC_summary_miRNAs_remout.txt", na = "NA", quote=FALSE, sep=",")

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