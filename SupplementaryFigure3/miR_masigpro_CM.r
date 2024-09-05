#!/usr/bin/env Rscript
#========================== package installation scripts ==============================
if (!require("maSigPro")){
	BiocManager::install("maSigPro")
}
if (!require("UpSetR")){
	install.packages("UpSetR")
}
library("UpSetR")
library("maSigPro")
library("pheatmap")
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)


setwd("/projects/amitimeseries/work/smallRNA_AMI/miRBase_counts/")
paste("Current working directory:", getwd())
#=================== Masigpro call for differential analysis ================================


# this function  (i) saves the table for DE miRNAs
# (ii) plots the intersection of the mature and precursor
# miRNAs (detected, significant, DE)


compute_sig <- function(sample_file, cell_type, norm_exp_file){

  # sample_file = "HC/samples_HC_new.txt"
  # cell_type = "CM"
  # norm_exp_file = "EC/EC_mature_normalized_CPM.txt"    #"FB/hairpin_normalized_CPM.txt"

  samples_FB <- read.table(file.path(sample_file), header=TRUE)

  rownames(samples_FB) = c(paste("ctrl", rep(cell_type,4), 1:4, sep="_"),
                           paste("d1", rep(cell_type,4), 1:4, sep="_"),
                           paste("d3", rep(cell_type,4), 1:4, sep="_"),
                           paste("d7", rep(cell_type,4), 1:4, sep="_"),
                           paste("d14", rep(cell_type,4), 1:4, sep="_"),
                           paste("d28", rep(cell_type,4), 1:4, sep="_"))    # head(samples_FB)

  data.mir.FB = read.table(file.path(norm_exp_file), header=TRUE)  # head(data.mir.FB,5)

  design <- make.design.matrix(samples_FB, degree = 5)
  colnames(data.mir.FB) = rownames(samples_FB)
  data.mir.FBn = data.mir.FB[rowSums(data.mir.FB) > 25,]

  fitFB <- p.vector(data = data.mir.FBn, design, Q = 0.01, MT.adjust = "BH", min.obs = 10)

  tFB <- T.fit(data = fitFB, step.method = "forward", alfa=0.01)

  get <- get.siggenes(tFB, rsq = 0.6, vars="groups")

  # see.genes(get$sig.genes$Group.1vsGroup.0, dis =design$dis,
  #           cluster.method="hclust" ,cluster.data = 1, k.mclust = TRUE)
  return(c(fitFB, get))
}

#---------------------------------------- plot heatmap -------------------------------------
# dev.off()
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

plot_heatmap <- function(data2plot){
  # data2plot = fitFB$SELEC
  data_subset_norm <- t(apply(data2plot, 1, cal_z_score))
  pheatmap(data_subset_norm,
           scale="row",
           clustering_distance_rows = "correlation", cluster_rows = TRUE,cutree_rows=5,
           cluster_cols=F,
           show_rownames=T,
           border_color=NA,
           filename="Heatmap_significant_miR_CM_mod.pdf"
  )
}

#---------------------------------------- plot upset --------------------------------------
plot_upset <- function(obj){

  listInput <- list(d1vsd0 = obj$summary$Group.1vsGroup.0,
                    d3vsd0 = obj$summary$Group.3vsGroup.0,
                    d7vsd0 = obj$summary$Group.7vsGroup.0,
                    d14vsd0 = obj$summary$Group.14vsGroup.0,
                    d28vsd0 = obj$summary$Group.28vsGroup.0)

  upsobj = UpSetR::upset(fromList(listInput), order.by = "freq",
                        sets.x.label = "DE miRNAs per timepoint",
                        text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 0))
  return(upsobj)
}


# # ==================================== Cardiomyocytes ======================================
cm = compute_sig(sample_file = "./configurations/samples_HC_new.txt",
                 cell_type = "CM",
                 norm_exp_file = "./expression/CM/CM_mature_normalized_CPM.txt")

# upset plot
pdf(file="upset_CM_miR.1.pdf", width = 4, height = 6, onefile=F) # or other device
plot_upset(cm)

dev.off()
