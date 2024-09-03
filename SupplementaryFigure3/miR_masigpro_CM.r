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
# library(ComplexUpset)
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

  sample_file = "HC/samples_HC_new.txt"
  cell_type = "CM"
  norm_exp_file = "EC/EC_mature_normalized_CPM.txt"    #"FB/hairpin_normalized_CPM.txt"

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
           border_color=NA# filename="Heatmap_significant_miR_FB_mod.png"
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
#---------------------------- pool expression matrix for stats analysis -----------------------
gen_pre_3p_5p_expr <- function(premature_5p_3p, ind3p, ind5p, ind_pre, tmp){
  # rm(df3p, df5p , dfpre)   ind3p = 52; ind5p = 28; ind_pre = 3; tmp = "d0"
  df3p = premature_5p_3p %>% select(mir3p, c(ind3p:(ind3p+3)))     #ctrl_FB_1, ctrl_FB_2, ctrl_FB_3, ctrl_FB_4)
  df5p = premature_5p_3p %>% select(mir5p, c(ind5p:(ind5p+3)))     #ctrl_FB_1.y, ctrl_FB_2.y, ctrl_FB_3.y, ctrl_FB_4.y)
  dfpre = premature_5p_3p %>% select(premature, c(ind_pre:(ind_pre+3)))         #ctrl_FB_1.x, ctrl_FB_2.x, ctrl_FB_3.x, ctrl_FB_4.x)
  
  return(assign(paste0(tmp,"_pre_3p_5p"),cbind(melt(df3p, value.name = paste0("three_expr")),
                                               melt(dfpre, value.name = paste0("pre_expr")),
                                               melt(df5p, value.name = paste0("five_expr")), timepoint = tmp)))
}
#------------------------------------ plot kruskal stats ------------------------------------
library(rstatix)
plot_stats <- function(premature_5p_3p, celltype, comparison, method, metric){
  #
  # 'comparison' is the text to appear on the plots
  # 'method' is the method for FD correction
  # 'metric' is the score to be tested for
  # premature_5p_3p=fb_premature_5p_3p; celltype = "FB"; comparison="all DE miRNAs"; method="none"; metric="log_ratio_prev"

  d0_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
  d1_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
  d3_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
  d7_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
  d14_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
  d28_pre_3p_5p = gen_pre_3p_5p_expr(premature_5p_3p, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
  pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p ,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)
  
  #Metrics
  # pre_3p_5p$loglog_ratio = log(log(abs(pre_3p_5p$three_expr - pre_3p_5p$five_expr))) #* (pre_3p_5p$pre_expr))
  pre_3p_5p$log_ratio = log(pre_3p_5p$three_expr) - log(pre_3p_5p$five_expr) #* (pre_3p_5p$pre_expr))
  pre_3p_5p$log_abs_diff = log(abs(pre_3p_5p$three_expr - pre_3p_5p$five_expr)) #* (pre_3p_5p$pre_expr))
  # pre_3p_5p$log_maxratio = log(max(c(pre_3p_5p$three_expr,pre_3p_5p$five_expr)) / pre_3p_5p$three_expr + pre_3p_5p$five_expr) #* (pre_3p_5p$pre_expr))
  pre_3p_5p$ratio = (pre_3p_5p$three_expr /  pre_3p_5p$five_expr)
  pre_3p_5p_sel = pre_3p_5p[,c("timepoint", metric, "log_ratio")] #as.data.frame(pre_3p_5p[, c(10,11)] ) # # #head(pre_3p_5p)
  pre_3p_5p_sel <- pre_3p_5p_sel %>% reorder_levels(timepoint, order = c("d0", "d1", "d3", "d7", "d14", "d28"))
  colnames(pre_3p_5p_sel) = c("timepoint", "metric", "log_ratio")
  
  log_ratio_plot <- ggboxplot(pre_3p_5p_sel, x = "timepoint", 
                              y = "log_ratio", fill = "timepoint", 
                              ylab="log_ratio",  ylim = c(-11, 6)) + 
  labs(
    title = paste0(celltype,": ","within DE miRNAs for ",comparison)
  ) + theme(legend.position = "none")

  res.kruskal <- pre_3p_5p_sel %>% kruskal_test(metric ~ timepoint)
  
  # Pairwise comparisons using Dunn’s test:
  stat.test <- pre_3p_5p_sel %>% dunn_test(metric ~ timepoint,
                                           p.adjust.method = method)#"bonferroni", BH, "none")
  
  # Visualization: box plots with stats n p-values
  pre_3p_5p_sel$timepoint = as.factor(pre_3p_5p_sel$timepoint)
  stat.test <- stat.test %>% add_xy_position(x = "timepoint") # add_y_position()
  bxp <- ggboxplot(pre_3p_5p_sel, x = "timepoint", y = "metric", fill = "timepoint", ylab=metric) + 
    stat_pvalue_manual(stat.test, hide.ns = TRUE) +
    labs(
      title = paste0(celltype,": ","within DE miRNAs for ",comparison),
      subtitle = get_test_label(res.kruskal, detailed = TRUE),
      caption = get_pwc_label(stat.test, type = "expression")
    ) + theme(legend.position = "none")
  
  figure <- ggarrange(log_ratio_plot,bxp,
                      labels = c("i", "ii"),
                      ncol = 1, nrow = 2)
  print(figure)
}
# ------------------------------- scatterplot for EDA ------------------------------

plot_scatter <- function(pre_3p_5p, x_axis, bool_conf, leg_start){
  # Explorative scatter plot of 3p and 5p expression correlation
  scp_data <- (pre_3p_5p %>% select(three_expr, five_expr, timepoint)
               %>% filter(timepoint == "d0"|
                            timepoint == "d1"|
                            timepoint == "d3"|
                            timepoint == "d14"|
                            timepoint == "d28"))
  scp_data <- scp_data %>%
    reorder_levels(timepoint, order = c("d0", "d1", "d3", "d7", "d14", "d28"))
  
  #https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
  scp_data$logthree_expr <- log(scp_data$three_expr)
  scp_data$logfive_expr <- log(scp_data$five_expr)
  scp_data$logratio <- log(scp_data$three_expr) - log(scp_data$five_expr)
  sp <- ggscatter(scp_data, x = x_axis, y = "logratio",
                  color = "timepoint", palette = "jco", size = 1,
                  add = "reg.line",
                  # add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = bool_conf)
  sp + stat_cor(method = "spearman",
                aes(color = timepoint), label.x = leg_start, size = 3.5)
  
}

# ------------------------------- to merge the data frames --------------------------------
build_df <- function(df1, df2, miR2probe){
  
  # miR2probe = ec$sig.genes$Group.14vsGroup.0$sig.profiles
  
  mature = as.data.frame(miR2probe)
  premature = str_replace(rownames(miR2probe), "-5p", "") %>% str_replace("-3p", "")
  
  premature_5p = merge(as.data.frame(premature), df1, by.x = "premature", by.y = 0)     #View(premature_5p)
  
  premature_5p$mir5p = paste0(premature_5p$premature, "-5p")
  premature_5p$mir3p = paste0(premature_5p$premature, "-3p")
  
  premature_5p_3p = merge(premature_5p, as.data.frame(df2), by.x="mir3p", by.y=0) %>%
    merge(as.data.frame(df2), by.x="mir5p", by.y=0)
  return(premature_5p_3p)
}
#====================================== Fibroblasts ====================================== 
# fb_pre = compute_sig(sample_file = "HC/samples_HC_new.txt",
#                  cell_type = "FB",
#                  norm_exp_file = "FB/hairpin_normalized_CPM.txt")

# # get the significant premature miRNA id and plot the heatmap for their 5p and 3p forms

# fb = compute_sig(sample_file = "HC/samples_HC_new.txt",
#                  cell_type = "FB",
#                  norm_exp_file = "FB/FB_mature_normalized_CPM.txt")

# # plot_heatmap(fb$sig.genes$Group.7vsGroup.0$sig.profiles)

# # build a function to return the merged larger dataframe with the premature expression
# # and corresponding 5p and 3p mature
# file = "FB/FB_mature_normalized_CPM.txt"
# df2 = read.table(file.path(file), header=TRUE) #fb  #the mature file

# file = "FB/hairpin_normalized_CPM.txt"
# df1 = read.table(file.path(file), header=TRUE) #fb  #the mature file

# #Case 3: checking usage of DE premature miRNAs at each timepoint
# rm(fb_premature_5p_3p, a, b, c, d, e, f)

# de_mir = rbind(fb$sig.genes$Group.1vsGroup.0$sig.profiles,
#                fb$sig.genes$Group.3vsGroup.0$sig.profiles,
#                fb$sig.genes$Group.7vsGroup.0$sig.profiles,
#                fb$sig.genes$Group.14vsGroup.0$sig.profiles,
#                fb$sig.genes$Group.28vsGroup.0$sig.profiles)
# fb_premature_5p_3p = build_df(df1, df2, de_mir) #<--
# a = plot_stats(fb_premature_5p_3p, "FB","all DE miRNAs", "BH", "log_abs_diff")
# # case 4
# rm(fb_premature_5p_3p)
# fb_premature_5p_3p = build_df(df1, df2,fb$sig.genes$Group.1vsGroup.0$sig.profiles) #
# b = plot_stats(fb_premature_5p_3p, "FB", "D1 vs D0", "BH", "log_abs_diff")
# fb_premature_5p_3p = build_df(df1, df2,fb$sig.genes$Group.3vsGroup.0$sig.profiles) #
# c = plot_stats(fb_premature_5p_3p, "FB", "D3 vs D0", "BH", "log_abs_diff")
# fb_premature_5p_3p = build_df(df1, df2,fb$sig.genes$Group.7vsGroup.0$sig.profiles)
# d = plot_stats(fb_premature_5p_3p,  "FB","D7 vs D0", "BH", "log_abs_diff")
# fb_premature_5p_3p = build_df(df1, df2,fb$sig.genes$Group.14vsGroup.0$sig.profiles)#<--
# e = plot_stats(fb_premature_5p_3p,  "FB","D14 vs D0", "BH", "log_abs_diff")
# fb_premature_5p_3p = build_df(df1, df2,fb$sig.genes$Group.28vsGroup.0$sig.profiles)
# f = plot_stats(fb_premature_5p_3p, "FB", "D28 vs D0", "BH", "log_abs_diff")

# figure <- ggarrange(a,b,c, #d,e,f, #
#                     labels = c("A", "B", "C"), #, "D","E","F"
#                     ncol = 3, nrow = 1)
# figure

# # .x: premature, .y: 5p; .null: 3p
# # perform correlation of the 3p and pri expression
# library(dplyr); library(reshape2) ;



# dev.off()
# fb_premature_5p_3p = build_df(fb_pre$SELEC, fb$SELEC, fb$SELEC)

# d0_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
# d1_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
# d3_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
# d7_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
# d14_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
# d28_pre_3p_5p = gen_pre_3p_5p_expr(fb_premature_5p_3p, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
# pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p ,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)

# w=plot_scatter(pre_3p_5p,"logthree_expr", bool_conf=TRUE, leg_start=1)
# x=plot_scatter(pre_3p_5p,"logfive_expr", TRUE, leg_start=6)
# y=plot_scatter(pre_3p_5p,"logthree_expr", FALSE, leg_start=1)
# z=plot_scatter(pre_3p_5p,"logfive_expr", FALSE, leg_start=6)

# figure <- ggarrange(w,x,y,z,
#                     labels = c("A", "B", "C", "D"),
#                     ncol = 2, nrow = 2)
# figure


# #----------------------------draw boxplots for EDA--------------------------------------
# drawboxplot <- function(d0_pre_3p_5p, tmp){

#   boxplot(log(d0_pre_3p_5p$pre_expr),
#           log(d0_pre_3p_5p$three_expr),
#           log(d0_pre_3p_5p$five_expr),
#           log(abs(d0_pre_3p_5p$three_expr - d0_pre_3p_5p$five_expr)),
#           log(d0_pre_3p_5p$three_expr) - log(d0_pre_3p_5p$five_expr),
#           names = c("pre_expr", "3p", "5p", "3p - 5p", "3p/5p"),
#           col = c(2,3,4,5,6),
#           main=paste0("Fibroblasts:",tmp))
#       # , ylim=c(-13,11)
# }

# drawboxplot(d0_pre_3p_5p, "D0")
# drawboxplot(d1_pre_3p_5p, "D1")
# drawboxplot(d3_pre_3p_5p, "D3")
# drawboxplot(d7_pre_3p_5p, "D7")
# drawboxplot(d14_pre_3p_5p, "D14")
# drawboxplot(d28_pre_3p_5p, "D28")
# # log of ratio
# boxplot(log(d0_pre_3p_5p$three_expr) - log(d0_pre_3p_5p$five_expr),
#         log(d1_pre_3p_5p$three_expr) - log(d1_pre_3p_5p$five_expr),
#         log(d3_pre_3p_5p$three_expr) - log(d3_pre_3p_5p$five_expr),
#         log(d7_pre_3p_5p$three_expr) - log(d7_pre_3p_5p$five_expr),
#         log(d14_pre_3p_5p$three_expr) - log(d14_pre_3p_5p$five_expr),
#         log(d28_pre_3p_5p$three_expr) - log(d28_pre_3p_5p$five_expr))
# # log of difference over premature
# boxplot(
#   log10(abs(d0_pre_3p_5p$three_expr - d0_pre_3p_5p$five_expr) / abs(d0_pre_3p_5p$pre_expr)),
#   log10(abs(d1_pre_3p_5p$three_expr - d1_pre_3p_5p$five_expr) / abs(d1_pre_3p_5p$pre_expr)),
#   log10(abs(d3_pre_3p_5p$three_expr - d3_pre_3p_5p$five_expr)/ abs(d1_pre_3p_5p$pre_expr)),
#   log10(abs(d7_pre_3p_5p$three_expr - d7_pre_3p_5p$five_expr)/ abs(d1_pre_3p_5p$pre_expr)),
#   log10(abs(d14_pre_3p_5p$three_expr - d14_pre_3p_5p$five_expr)/ abs(d1_pre_3p_5p$pre_expr)),
#   log10(abs(d28_pre_3p_5p$three_expr - d28_pre_3p_5p$five_expr)/ abs(d1_pre_3p_5p$pre_expr)),
#   main = "log of difference over premature (significant mature miRNAs)",
#   names = c("D0", "D1", "D3", "D7", "D14", "D28"),
#   col = c(2,3,4,5,6,7), ylim=c(-2,2)
# )

# # log of difference
# boxplot(
#   log10(abs(d0_pre_3p_5p$three_expr - d0_pre_3p_5p$five_expr) ),
#   log10(abs(d1_pre_3p_5p$three_expr - d1_pre_3p_5p$five_expr) ),
#   log10(abs(d3_pre_3p_5p$three_expr - d3_pre_3p_5p$five_expr)),
#   log10(abs(d7_pre_3p_5p$three_expr - d7_pre_3p_5p$five_expr)),
#   log10(abs(d14_pre_3p_5p$three_expr - d14_pre_3p_5p$five_expr)),
#   log10(abs(d28_pre_3p_5p$three_expr - d28_pre_3p_5p$five_expr)),
#   main = "log of difference (significant mature miRNAs)",
#   names = c("D0", "D1", "D3", "D7", "D14", "D28"),
#   col = c(2,3,4,5,6,7), ylim=c(-2,5)
# )


# #Alternative way to compute arm usage ratio
# # ctrlratio = fb_premature_5p_3p$ctrl_FB_1 / (fb_premature_5p_3p$ctrl_FB_1.y + fb_premature_5p_3p$ctrl_FB_1)
# # d1ratio = fb_premature_5p_3p$d1_FB_1 / (fb_premature_5p_3p$d1_FB_1.y + fb_premature_5p_3p$d1_FB_1)
# # d3ratio = fb_premature_5p_3p$d3_FB_1 / (fb_premature_5p_3p$d3_FB_1.y + fb_premature_5p_3p$d3_FB_1)
# # d7ratio = fb_premature_5p_3p$d7_FB_1 / (fb_premature_5p_3p$d7_FB_1.y + fb_premature_5p_3p$d7_FB_1)
# # d14ratio = fb_premature_5p_3p$d14_FB_1 / (fb_premature_5p_3p$d14_FB_1.y + fb_premature_5p_3p$d14_FB_1)
# # d28ratio = fb_premature_5p_3p$d28_FB_1 / (fb_premature_5p_3p$d28_FB_1.y + fb_premature_5p_3p$d28_FB_1)
# # .x: premature, .y: 5p; .null: 3p
# ctrlratio = fb_premature_5p_3p$ctrl_FB_1 / (fb_premature_5p_3p$ctrl_FB_1.y)  # + fb_premature_5p_3p$ctrl_FB_1)
# d1ratio = fb_premature_5p_3p$d1_FB_1 / (fb_premature_5p_3p$d1_FB_1.y)  # + fb_premature_5p_3p$d1_FB_1)
# d3ratio = fb_premature_5p_3p$d3_FB_1 / (fb_premature_5p_3p$d3_FB_1.y)  # + fb_premature_5p_3p$d3_FB_1)
# d7ratio = fb_premature_5p_3p$d7_FB_1 / (fb_premature_5p_3p$d7_FB_1.y)  # + fb_premature_5p_3p$d7_FB_1)
# d14ratio = fb_premature_5p_3p$d14_FB_1 / (fb_premature_5p_3p$d14_FB_1.y)# + fb_premature_5p_3p$d14_FB_1)
# d28ratio = fb_premature_5p_3p$d28_FB_1 / (fb_premature_5p_3p$d28_FB_1.y)  # + fb_premature_5p_3p$d28_FB_1)

# pred1ratio = fb_premature_5p_3p$d1_FB_1.x / fb_premature_5p_3p$ctrl_FB_1.x
# pred3ratio = fb_premature_5p_3p$d3_FB_1.x / fb_premature_5p_3p$ctrl_FB_1.x
# pred7ratio = fb_premature_5p_3p$d7_FB_1.x / fb_premature_5p_3p$ctrl_FB_1.x
# pred14ratio = fb_premature_5p_3p$d14_FB_1.x / fb_premature_5p_3p$ctrl_FB_1.x
# pred28ratio = fb_premature_5p_3p$d28_FB_1.x / fb_premature_5p_3p$ctrl_FB_1.x

# ratio_results = data.frame(premiRNAs = fb_premature_5p_3p$pre,
#                             ctrl_usage=ctrlratio,
#                             d1_usage = d1ratio,
#                             d3_usage = d3ratio,
#                             d7_usage = d7ratio,
#                            d14_usage = d14ratio,
#                            d28_usage = d28ratio,
#                            d1_switch = (d1ratio / ctrlratio) / (pred1ratio),
#                            d3_switch = (d3ratio / ctrlratio)/(pred3ratio),
#                            d7_switch = (d7ratio / ctrlratio)/(pred7ratio),
#                            d14_switch = (d14ratio / ctrlratio) / (pred14ratio),
#                            d28_switch = (d28ratio / ctrlratio) / (pred28ratio))

# View(ratio_results)


# data_subset <- fb$SELEC[rownames(fb$SELEC) %in% select_pre, ]#head(data_subset)
# plot_heatmap(data_subset)

# #---------------------- plot and save upset plot (FB) ----------------------
# png(file="upset_FB_premature_miR_v1.png") # or other device
# plot_upset(fb_pre)
# grid.text("DE premarture miRNAs in Fibroblasts", x = 0.65, y = 0.95,
#                     gp = gpar(fontsize = 15))
# dev.off()
# # ==================================== Cardiomyocytes ======================================
cm = compute_sig(sample_file = "HC/samples_HC_new.txt",
                 cell_type = "CM",
                 norm_exp_file = "CM/CM_mature_normalized_CPM.txt")

# # cm_pre = compute_sig(sample_file = "HC/samples_HC_new.txt",
# #                  cell_type = "CM",
# #                  norm_exp_file = "CM/hairpin_normalized_CPM.txt")

# upset plot
pdf(file="upset_CM_miR.1.pdf", width = 4, height = 6, onefile=F) # or other device
plot_upset(cm)
# grid.text("DE miRNAs in Cardiomyocytes", x = 0.65, y = 0.95,
#           gp = gpar(fontsize = 15))
dev.off()

# # ---------------------------- build the dataframe -----------------------------
# # build a function to return the merged larger dataframe with the premature expression
# # and corresponding 5p and 3p mature
# file = "CM/CM_mature_normalized_CPM.txt"
# df2 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# file = "CM/hairpin_normalized_CPM.txt"
# df1 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# #Case 3: checking usage of DE premature miRNAs at each timepoint
# rm(cm_premature_5p_3p, a,b,c,d,e,f)

# de_mir = rbind(cm$sig.genes$Group.1vsGroup.0$sig.profiles,
#       cm$sig.genes$Group.3vsGroup.0$sig.profiles,
#       cm$sig.genes$Group.7vsGroup.0$sig.profiles,
#       cm$sig.genes$Group.14vsGroup.0$sig.profiles,
#       cm$sig.genes$Group.28vsGroup.0$sig.profiles)
# cm_premature_5p_3p = build_df(df1, df2, de_mir); #<--
# a = plot_stats(cm_premature_5p_3p, "CM", "all DE miRNAs", "BH", "log_abs_diff")
# # case 4
# rm(cm_premature_5p_3p)
# cm_premature_5p_3p = build_df(df1, df2,cm$sig.genes$Group.1vsGroup.0$sig.profiles) #
# b = plot_stats(cm_premature_5p_3p, "CM","D1 vs D0", "BH", "log_abs_diff")
# cm_premature_5p_3p = build_df(df1, df2, cm$sig.genes$Group.3vsGroup.0$sig.profiles) #
# c = plot_stats(cm_premature_5p_3p, "CM", "D3 vs D0", "BH", "log_abs_diff")
# cm_premature_5p_3p = build_df(df1, df2, cm$sig.genes$Group.7vsGroup.0$sig.profiles)
# d = plot_stats(cm_premature_5p_3p, "CM", "D7 vs D0", "BH", "log_abs_diff")
# cm_premature_5p_3p = build_df(df1, df2,cm$sig.genes$Group.14vsGroup.0$sig.profiles)#<--
# e = plot_stats(cm_premature_5p_3p, "CM", "D14 vs D0", "BH", "log_abs_diff")
# cm_premature_5p_3p = build_df(df1, df2,cm$sig.genes$Group.28vsGroup.0$sig.profiles)
# f = plot_stats(cm_premature_5p_3p, "CM", "D28 vs D0", "BH", "log_abs_diff")

# figure <- ggarrange( a,b,c,#d,e,f,#
#                     labels = c("A","B","C"),
#                     ncol = 3, nrow = 1)
# figure

# # ------------------------------------------- scatter plot -----------------------------------------
# dev.off()
# cm_premature_5p_3p = build_df(cm_pre$SELEC, cm$SELEC, cm$SELEC)

# d0_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
# d1_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
# d3_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
# d7_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
# d14_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
# d28_pre_3p_5p = gen_pre_3p_5p_expr(cm_premature_5p_3p, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
# pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p ,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)

# w=plot_scatter(pre_3p_5p,"logthree_expr", bool_conf=TRUE, leg_start=1)
# x=plot_scatter(pre_3p_5p,"logfive_expr", TRUE, leg_start=6)
# y=plot_scatter(pre_3p_5p,"logthree_expr", FALSE, leg_start=1)
# z=plot_scatter(pre_3p_5p,"logfive_expr", FALSE, leg_start=6)

# figure <- ggarrange(w,x,y,z,
#                     labels = c("A", "B", "C", "D"),
#                     ncol = 2, nrow = 2)
# figure
# # ======================================= Endothelial Cells ======================================
# ec = compute_sig(sample_file = "HC/samples_HC_new.txt",
#                  cell_type = "EC",
#                  norm_exp_file = "EC/EC_mature_normalized_CPM.txt")

# ec_pre = compute_sig(sample_file = "HC/samples_HC_new.txt",
#                  cell_type = "EC",
#                  norm_exp_file = "EC/hairpin_normalized_CPM.txt")

# png(file="upset_EC_premature_miR_v1.png") # or other device
# plot_upset(ec_pre)
# grid.text("DE premature miRNAs in Endothelial Cells", x = 0.65, y = 0.95,
#           gp = gpar(fontsize = 15))
# dev.off()

# # ---------------------------- build the dataframe for stats -----------------------------
# # build a function to return the merged larger dataframe with the premature expression
# # and corresponding 5p and 3p mature
# file = "EC/EC_mature_normalized_CPM.txt"
# df2 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# file = "EC/hairpin_normalized_CPM.txt"
# df1 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# #Case 3: checking usage of DE mature miRNAs in each timepoint
# rm(ec_premature_5p_3p, de_mir)
# de_mir = rbind(ec$sig.genes$Group.1vsGroup.0$sig.profiles,
#                ec$sig.genes$Group.3vsGroup.0$sig.profiles,
#                ec$sig.genes$Group.7vsGroup.0$sig.profiles,
#                ec$sig.genes$Group.14vsGroup.0$sig.profiles,
#                ec$sig.genes$Group.28vsGroup.0$sig.profiles)
# ec_premature_5p_3p = build_df(df1, df2,de_mir)
# a = plot_stats(ec_premature_5p_3p, "EC","all DE miRNAs", "BH", "log_abs_diff")

# #Case 4: checking usage of DE premature miRNAs at each timepoint
# rm(ec_premature_5p_3p)
# ec_premature_5p_3p = build_df(df1, df2, ec$sig.genes$Group.1vsGroup.0$sig.profiles)
# b=plot_stats(ec_premature_5p_3p, "EC","D1 vs D0","BH", "log_abs_diff")
# ec_premature_5p_3p = build_df(df1, df2, ec$sig.genes$Group.3vsGroup.0$sig.profiles) #*
# c=plot_stats(ec_premature_5p_3p,"EC", "D3 vs D0","BH", "log_abs_diff")
# ec_premature_5p_3p = build_df(df1, df2, ec$sig.genes$Group.7vsGroup.0$sig.profiles)
# d=plot_stats(ec_premature_5p_3p,"EC", "D7 vs D0", "none", "log_abs_diff")
# ec_premature_5p_3p = build_df(df1, df2, ec$sig.genes$Group.14vsGroup.0$sig.profiles)
# e = plot_stats(ec_premature_5p_3p,"EC", "D14 vs D0", "none", "log_abs_diff")
# ec_premature_5p_3p = build_df(df1, df2, ec$sig.genes$Group.28vsGroup.0$sig.profiles)
# f=plot_stats(ec_premature_5p_3p, "EC","D28 vs D0", "none", "log_abs_diff")

# figure <- ggarrange(a,b,c, #d,e,f,#
#                     labels = c("A", "B", "C"),# "D","E","F", "G","H", "I"
#                     ncol = 3, nrow = 1)
# figure
# # ---------------------------build dataframe for scatterplot -------------------------------
# rm(pre_3p_5p, pre_3p_5p_sel)
# cm_premature_5p_3p = build_df(ec_pre$SELEC, ec$SELEC, ec$SELEC)
# d0_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
# d1_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
# d3_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
# d7_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
# d14_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
# d28_pre_3p_5p = gen_pre_3p_5p_expr(ec_premature_5p_3p, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
# pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p ,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)

# # --------------------------------- scatter plot -----------------------------------
# dev.off()
# w=plot_scatter(pre_3p_5p,"logthree_expr", bool_conf=TRUE, leg_start = -2)
# x=plot_scatter(pre_3p_5p,"logfive_expr", bool_conf=TRUE, leg_start = 6)
# y=plot_scatter(pre_3p_5p,"logthree_expr", bool_conf=FALSE, leg_start = -2)
# z=plot_scatter(pre_3p_5p,"logfive_expr", bool_conf=FALSE, leg_start = 6)

# figure <- ggarrange(w,x,y,z,
#                     labels = c("A", "B", "C", "D"),
#                     ncol = 2, nrow = 2)
# figure
# # ======================================== Hematopoietic cells =========================================

# hc = compute_sig(sample_file = "HC/samples_HC_new.txt",
#                  cell_type = "HC",
#                  norm_exp_file = "HC/HC_mature_normalized_CPM.txt")

# # hc_pre = compute_sig(sample_file = "HC/samples_HC_new.txt",
# #                  cell_type = "HC",
# #                  norm_exp_file = "HC/hairpin_normalized_CPM.txt")

# # png(file="upset_HC_premature_miR_v1.png") # or other device
# # plot_upset(hc_pre)
# # grid.text("DE premature miRNAs in Hematopoietic Cells", x = 0.65, y = 0.95,
# #           gp = gpar(fontsize = 15))
# # dev.off()

# # ------------------------------------------- build the dataframe -----------------------------------------
# file = "HC/HC_mature_normalized_CPM.txt"
# df2 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# file = "HC/hairpin_normalized_CPM.txt"
# df1 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# #Case 3: checking usage of all DE mature miRNAs in each timepoint
# rm(hc_premature_5p_3p, a, b, c, d, e, f)
# de_mir = rbind(hc$sig.genes$Group.1vsGroup.0$sig.profiles,
#                hc$sig.genes$Group.3vsGroup.0$sig.profiles,
#                hc$sig.genes$Group.7vsGroup.0$sig.profiles,
#                hc$sig.genes$Group.14vsGroup.0$sig.profiles,
#                hc$sig.genes$Group.28vsGroup.0$sig.profiles)
# hc_premature_5p_3p = build_df(df1, df2,de_mir)
# a = plot_stats(hc_premature_5p_3p, "HC", "all DE miRNAs", "BH", "log_abs_diff") #

# #Case 4: checking usage of DE premature miRNAs at each timepoint
# rm(ec_premature_5p_3p)
# hc_premature_5p_3p = build_df(df1, df2, hc$sig.genes$Group.1vsGroup.0$sig.profiles)
# b=plot_stats(hc_premature_5p_3p, "HC","D1 vs D0", "BH", "log_abs_diff")
# hc_premature_5p_3p = build_df(df1, df2, hc$sig.genes$Group.3vsGroup.0$sig.profiles) #*
# c=plot_stats(hc_premature_5p_3p,"HC", "D3 vs D0", "BH", "log_abs_diff")
# hc_premature_5p_3p = build_df(df1, df2, hc$sig.genes$Group.7vsGroup.0$sig.profiles)
# d=plot_stats(hc_premature_5p_3p, "HC","D7 vs D0", "BH", "log_abs_diff")
# hc_premature_5p_3p = build_df(df1, df2, hc$sig.genes$Group.14vsGroup.0$sig.profiles)
# e = plot_stats(hc_premature_5p_3p, "HC", "D14 vs D0", "BH", "log_abs_diff")
# hc_premature_5p_3p = build_df(df1, df2, hc$sig.genes$Group.28vsGroup.0$sig.profiles)
# f=plot_stats(hc_premature_5p_3p,  "HC", "D28 vs D0",  "BH", "log_abs_diff")

# figure <- ggarrange(a,b,c,#d,e,f, # 
#                     labels = c("A","B","C", "G","H", "I"),
#                     ncol = 3, nrow = 1)
# figure
# # ---------------------------------------------- scatter plot -------------------------------------------
# dev.off()
# hc_premature_5p_3p = build_df(hc_pre$SELEC, hc$SELEC, hc$SELEC)

# d0_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
# d1_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
# d3_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
# d7_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
# d14_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
# d28_pre_3p_5p = gen_pre_3p_5p_expr(hc_premature_5p_3p, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
# pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p ,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)

# w=plot_scatter(pre_3p_5p,"logthree_expr", bool_conf=TRUE, leg_start=1)
# x=plot_scatter(pre_3p_5p,"logfive_expr", TRUE, leg_start=6)
# y=plot_scatter(pre_3p_5p,"logthree_expr", FALSE, leg_start=1)
# z=plot_scatter(pre_3p_5p,"logfive_expr", FALSE, leg_start=6)

# figure <- ggarrange(w,x,y,z,
#                     labels = c("A", "B", "C", "D"),
#                     ncol = 2, nrow = 2)
# figure
# # ---------------------------- stats per miRNA ------------------------------------------

# # perform mann whitney test for each significant miRNA
# perform_mann <- function(df1, df2, sig){
  
#   # sig=cm$SELEC
#   # Input: sig: dataframe with the expression matrix (of matched 3p and 5p)
#   df_exp = build_df(df1, df2, miR2probe = sig)
  
#   d0_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 28, ind5p = 52, ind_pre = 4, tmp = "d0")
#   d1_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 32, ind5p = 56, ind_pre = 8, tmp = "d1")
#   d3_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 36, ind5p = 60, ind_pre = 12, tmp = "d3")
#   d7_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 40, ind5p = 64, ind_pre = 16, tmp = "d7")
#   d14_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 44, ind5p = 68, ind_pre = 20, tmp = "d14")
#   d28_pre_3p_5p = gen_pre_3p_5p_expr(df_exp, ind3p = 48, ind5p = 72, ind_pre = 24, tmp = "d28")
#   pre_3p_5p = rbind(d0_pre_3p_5p, d1_pre_3p_5p,d3_pre_3p_5p, d7_pre_3p_5p, d14_pre_3p_5p, d28_pre_3p_5p)

#   #metrics to compare
#   pre_3p_5p$log_ratio = log(pre_3p_5p$three_expr) - log(pre_3p_5p$five_expr)
#   pre_3p_5p$log_abs_diff = log(abs(pre_3p_5p$three_expr - pre_3p_5p$five_expr))
#   pre_3p_5p$weight = pre_3p_5p$log_abs_diff
  
#   premature_list = as.character(unique(pre_3p_5p$premature)) # list of premature miRNAs
#   pre_3p_5p$timepoint <- ordered(pre_3p_5p$timepoint,
#                                     levels = c("d0", "d1", "d3", "d7", "d14", "d28"))
#   results = data.frame(); results_all = data.frame()
  
#   for (var in premature_list){
#     tmp1 = data.frame(premature = c(var))
#     results_all = rbind(results_all, tmp1)
#   }
  
#   for (lv in 2:length(levels(pre_3p_5p$timepoint))){
    
#     cmp_lvl = levels(pre_3p_5p$timepoint)[lv]
#     # cmp_lvl = "d14"
#     for (var in premature_list){
 
#       # var = "mmu-miR-10b"
#       df_exp_sel = subset(pre_3p_5p, premature == var & (timepoint == cmp_lvl | timepoint == "d0"))
#       df_exp_sel1 <- df_exp_sel %>% select(weight, timepoint)
#       df_exp_sel1$timepoint <- ordered(df_exp_sel1$timepoint, levels = c("d0", cmp_lvl))
#       rownames(df_exp_sel1) <- 1:nrow(df_exp_sel1)
      
#       sts = group_by(df_exp_sel, timepoint) %>%
#         summarise(
#           #count = n(), mean = mean(weight, na.rm = TRUE), sd = sd(weight, na.rm = TRUE),
#           median = median(log_ratio, na.rm = TRUE),
#           # IQR = IQR(weight, na.rm = TRUE)
#         )

#       #computation: add man whitney test p-value
#       miR.stat.test <- df_exp_sel1  %>%
#         wilcox_test(weight ~ timepoint, paired = TRUE) %>%
#         add_significance() #stat.test
      
#       cat(var, "\t",miR.stat.test$p,"\t",miR.stat.test$p.signif,"\n")
#       tmp = data.frame(
#         # premature = c(var),
#                        ini_ratio = sts$median[1],
#                        to_ratio = sts$median[2],
#                       sig = c(miR.stat.test$p), 
#                        status = c(miR.stat.test$p.signif))
      
#       names(tmp) = c(
#         # paste0(cmp_lvl, "pre"),
#                      paste0(cmp_lvl, "ini_ratio"),
#                      paste0(cmp_lvl, "to_ratio"),
#                      paste0(cmp_lvl, "sig"), 
#                      paste0(cmp_lvl, "status"))
#       results = rbind(results, tmp)
      
#     } #end var loop
#     results$padj = p.adjust(results[,3], method = "BH")
#     names(results)[5] = paste0(cmp_lvl, "padj")
#     results_all = cbind(results_all,results)
   
#     results = data.frame()  
#   } #end lvl loop
  
#   return(results_all)
# }

# file = "CM/CM_mature_normalized_CPM.txt"
# df2 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# file = "CM/hairpin_normalized_CPM.txt"
# df1 = read.table(file.path(file), header=TRUE) #cm  #the mature file

# # CM_check_out = perform_mann(df1, df2, sig=cm$SELEC)
# # write.table(CM_check_out, file = "sig_CM_arm_usage.txt", sep = "\t", quote = FALSE)
# # --------------------------------------------------------------------------------------------------
# # rm(to_plot)
# # to_plot <- CM_check_out %>% filter_at(vars(ends_with("status")), any_vars(. == "*")) %>%
# #   select(premature, d1ini_ratio,ends_with("to_ratio"))
# # rownames(to_plot) = to_plot$premature

# # to_plot$d1ini_ratio[to_plot$d1ini_ratio <= -10] <- -10
# # to_plot$logFC[df$pval >= 0.05] <- 0

# # plot_heatmap(to_plot)
# #==================== need to remove the outlier samples and run again d0, rep 2 ==================
# print("hello: Running miRNA differential expression analysis removing outlier FB samples")
# samples_FB <- read.table(file.path("FB/samples_FB_remout.txt"), header=TRUE)
# rownames(samples_FB) = c(paste("ctrl", rep("FB",1), 1, sep="_"),
#                          paste("ctrl", rep("FB",2), 3:4, sep="_"),
#                          paste("d1", rep("FB",4), 1:4, sep="_"),
#                          paste("d3", rep("FB",4), 1:4, sep="_"),
#                          paste("d7", rep("FB",4), 1:4, sep="_"),
#                          paste("d14", rep("FB",4), 1:4, sep="_"),
#                          paste("d28", rep("FB",4), 1:4, sep="_"))

# head(samples_FB,5)

# data.mir.FB = read.table(file.path("FB/FB_mature_normalized_CPM.txt"), header=TRUE)
# drops <- c("ctrl_FB_2")
# data.mir.FB = data.mir.FB[ , !(colnames(data.mir.FB) %in% drops)]
# head(data.mir.FB,3)

# design <- make.design.matrix(samples_FB, degree = 5)
# colnames(data.mir.FB) = rownames(samples_FB)
# data.mir.FBn = data.mir.FB[rowSums(data.mir.FB) > 10,]

# dim(data.mir.FB)
# dim(data.mir.FBn)

# fitFB <- p.vector(data = data.mir.FBn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

# paste("# of sig. miRs", fitFB$i) # returns the number of significant genes
# paste("alpha", head(fitFB$alfa))# gives p-value at the Q false discovery control level
# head(fitFB$SELEC) # is a matrix with the significant genes and their expression values
# head(fitFB$p.vector) # vector containing the computed p-values
# head(fitFB$p.adjusted)

# pheatmap(fitFB$SELEC,
#          scale="row",
#          clustering_distance_rows = "correlation",
#          cluster_cols=F,
#          show_rownames=F,
#          border_color=NA,
#          filename="Heatmap_significant_miR_FB_remout.png"
# )

# tFB <- T.fit(data = fitFB, step.method = "two.ways.forward")

# get<-get.siggenes(tFB, vars="groups")
# write.table(get$summary, "FB_summary_miRNAs_remout.txt", na = "NA", quote=FALSE, sep=",")

# # names(get$sig.genes)
# # names(get$sig.genes$Group.1vsGroup.0)
# # View(get$sig.genes$Group.1vsGroup.0$sig.pvalues)
# # View(get$sig.genes$Group.1vsGroup.0$coefficients)
# # names(get$summary$Group.1vsGroup.0)

# png("upset_FB_miR_remout.png", width = 600, height = 600, res=100)
# listInput <- list(d1vsd0 = get$summary$Group.1vsGroup.0,
#                   d3vsd0 = get$summary$Group.3vsGroup.0,
#                   d7vsd0 = get$summary$Group.7vsGroup.0,
#                   d14vsd0 = get$summary$Group.14vsGroup.0,
#                   d28vsd0 = get$summary$Group.28vsGroup.0)
# upset(fromList(listInput), order.by = "freq")
# dev.off()
# #========================================================================================
