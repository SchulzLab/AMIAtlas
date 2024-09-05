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

#==================== need to remove the outlier samples and run again ==================
print("hello: Running miRNA differential expression analysis removing outlier HC samples")
samples_HC <- read.table(file.path("./configurations/samples_HC_remout.txt"), header=TRUE)
rownames(samples_HC) = c(paste("ctrl", rep("HC",2), 1:2, sep="_"), 
                         paste("ctrl", rep("HC",1), 4, sep="_"),
                         paste("d1", rep("HC",4), 1:4, sep="_"),
                         paste("d3", rep("HC",4), 1:4, sep="_"),
                         paste("d7", rep("HC",4), 1:4, sep="_"),
                         paste("d14", rep("HC",4), 1:4, sep="_"),
                         paste("d28", rep("HC",4), 1:4, sep="_"))

head(samples_HC,5) # removed ctrl_HC_3

data.mir.HC = read.table(file.path("./expression/HC/HC_mature_normalized_CPM.txt"), header=TRUE)
drops <- c("ctrl_HC_3")
data.mir.HC = data.mir.HC[ , !(colnames(data.mir.HC) %in% drops)]
head(data.mir.HC,3)

design <- make.design.matrix(samples_HC, degree = 5)
colnames(data.mir.HC) = rownames(samples_HC)
data.mir.HCn = data.mir.HC[rowSums(data.mir.HC) > 10,]

dim(data.mir.HC)
dim(data.mir.HCn)

fitHC <- p.vector(data = data.mir.HCn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

paste("# of sig. miRs", fitHC$i) # returns the number of significant genes
paste("alpha", head(fitHC$alfa))# gives p-value at the Q false discovery control level
head(fitHC$SELEC) # is a matrix with the significant genes and their expression values
head(fitHC$p.vector) # vector containing the computed p-values
head(fitHC$p.adjusted)

pheatmap(fitHC$SELEC,
         scale="row",
         clustering_distance_rows = "correlation",
         cluster_cols=F,
         show_rownames=F,
         border_color=NA,
         filename="Heatmap_significant_miR_HC_remout.png"
)

tHC <- T.fit(data = fitHC, step.method = "two.ways.forward")

get<-get.siggenes(tHC, vars="groups")

write.table(get$summary, "HC_summary_miRNAs_remout.txt", na = "NA", quote=FALSE, sep=",")

pdf("./results/upset_HC_miR_remout.1.pdf", 
        width = 4, 
        height = 6, onefile=F)
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