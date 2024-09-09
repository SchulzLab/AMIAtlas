setwd("/Volumes/Elements_1/smallRNA_MI/miRge3/")
getwd()

## ========================= load libraries ===================== 
library(ggrepel)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(matrixStats)
library(reshape2)
library(readr)
library("dplyr")    # Data manipulation
library(ggpubr)
# ------------------- create a dataset-------------------
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)

# ========================== create a dataset(toy) ====================
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3),
#             rep("sorgho1" , 3) , rep("poacee1" , 3) , rep("banana1" , 3) , rep("triticum1" , 3))
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 8)
# value <- abs(rnorm(24 , 0 , 15))
# data1 <- data.frame(specie,condition,value)
# 
# # Stacked + percent
# ggplot(data, aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="fill", stat="identity")
# 
# # grouped
# ggplot(data1, aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="dodge", stat="identity")
# 
# # Small multiple
# ggplot(data1, aes(fill=condition, y=value, x=specie)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_viridis(discrete = T) +
#   ggtitle("Studying 4 species..") +
#   theme_ipsum() +
#   xlab("")

# ---------------------------Example on toy dataset-----------------------
# create the toy dataset
timepoint <- c(rep("D0" , 7), rep("D1" , 7), rep("D3" , 7), rep("D7" , 7), rep("D14" , 7), rep("D28" , 7))
RNAtype <- rep(c("rRNA" , "Hairpin miRNA" , "snoRNA", "primary tRNA", "mature miRNA", "other ncRNA", "mature tRNA") , 6)
value <- abs(rnorm(42 , 0 , 15))
data2 <- data.frame(timepoint,RNAtype,value)
data2$timepoint <- factor(data2$timepoint , 
                          levels= c("D0",  "D1", "D3", "D7", "D14", "D28" ))

# Stack plot on toy dataset
ggplot(data2, aes(fill=RNAtype, y=value, x=timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), 
            position=position_stack(1.0), vjust=-0.5, size=2.5, 
            color = "black",fontface = "bold", check_overlap = T) +
  ggtitle("Hematopoetic Cells") +
  theme_ipsum() + xlab("Timepoints") + ylab("Percentage of total reads")

# Percent stacked Plot on toy dataset
ggplot(data2, aes(fill=RNAtype, y=value, x=timepoint)) + 
  geom_bar(position="fill", stat="identity") +
  # scale_fill_viridis(discrete = T) +
  geom_text(aes(label = format(round(value, 2), nsmall = 2)), 
            position=position_fill(1.0), vjust=-0.5, size=2.5, 
            color = "black",fontface = "bold", check_overlap = T) +
  ggtitle("Hematopoetic Cells") +
  theme_ipsum() +
  xlab("Timepoints")

# --------------- trying pie plot ---------------------------


## -------------------- end of example on Toy dataset -----------------
# To Do: Add remaining category
count.data <- data.frame(
  class = c("1st", "2nd", "3rd", "Crew"),
  n = c(325, 285, 706, 885),
  prop = c(14.8, 12.9, 32.1, 40.2)
)
count.data
# Add label position
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white")+
  scale_fill_manual(values = mycols) +
  theme_void()
## ============================= Real Dataset =============================== 



#To Do: Check for MDS plot for choosing the samples, getting rid of the outlier samples
# Outlier samples: ctrl_HC_3, ctrl_FB_2, d28_EC_3

snc_RNAdis_plot <- function(file, title, cell_type, select){
  # file = "CM.annotation.report.csv"; title = "Cardiomyocytes"; cell_type="CM"
  EC_annotation_report <- read_csv(file) # head(HC_annotation_report) #colnames(HC_annotation_report)
  #choose the columns
  # [1]"Sample name(s)" [2]"Total Input Reads" [3]"Trimmed Reads (all)"  [4]"Trimmed Reads (unique)" [5]"All miRNA Reads" [6]"Filtered miRNA Reads" [7]"Unique miRNAs"         
  # [8]"Hairpin miRNAs" [9]"mature tRNA Reads" [10]"primary tRNA Reads" [11]"snoRNA Reads" [12]"rRNA Reads" [13]"ncRNA others" [14]"mRNA Reads" [15] "Remaining Reads"
  EC <- EC_annotation_report[, c(1, 2, 6, 8, 9, 10, 11, 12, 13, 15)] #head(EC)
  View(EC)
  #generate the median of the samples table after remove outlier samples
  if(cell_type == "FB"){
    median_reads = data.frame(d0 = colMedians(as.matrix(EC[c(1,3,4),-1])), 
                              d1 = colMedians(as.matrix(EC[9:12,-1])),
                              d3 = colMedians(as.matrix(EC[17:20,-1])),
                              d7 = colMedians(as.matrix(EC[21:24,-1])),
                              d14 = colMedians(as.matrix(EC[5:8,-1])),
                              d28 = colMedians(as.matrix(EC[13:16,-1])))
  } else if(cell_type == "EC"){
    median_reads = data.frame(d0 = colMedians(as.matrix(EC[1:4,-1])), 
                              d1 = colMedians(as.matrix(EC[9:12,-1])),
                              d3 = colMedians(as.matrix(EC[17:20,-1])),
                              d7 = colMedians(as.matrix(EC[21:24,-1])),
                              d14 = colMedians(as.matrix(EC[5:8,-1])),
                              d28 = colMedians(as.matrix(EC[c(13,14,16),-1])))
  } else if(cell_type == "HC"){     # ctrl_HC_3
    median_reads = data.frame(d0 = colMedians(as.matrix(EC[c(1,2,4),-1])), 
                              d1 = colMedians(as.matrix(EC[9:12,-1])),
                              d3 = colMedians(as.matrix(EC[17:20,-1])),
                              d7 = colMedians(as.matrix(EC[21:24,-1])),
                              d14 = colMedians(as.matrix(EC[5:8,-1])),
                              d28 = colMedians(as.matrix(EC[13:16,-1])))
  }else{
    median_reads = data.frame(d0 = colMedians(as.matrix(EC[1:4,-1])), #Check: median(c(5140099, 16780543, 34440326,  4599833))
                            d1 = colMedians(as.matrix(EC[9:12,-1])),
                            d3 = colMedians(as.matrix(EC[17:20,-1])),
                            d7 = colMedians(as.matrix(EC[21:24,-1])),
                            d14 = colMedians(as.matrix(EC[5:8,-1])),
                            d28 = colMedians(as.matrix(EC[13:16,-1])))
  }
  rownames(median_reads) = c("Total Reads","miRNA Reads","Hairpin miRNAs",
                             "mature tRNA Reads","primary tRNA Reads", "snoRNA Reads",
                             "rRNA Reads","ncRNA others","Remaining Reads")
  View(median_reads)

  perc_reads = rbind(
                      ((median_reads[2,] / (median_reads[1,] - median_reads[7,] ))*100), #miRNA reads
                     ((median_reads[3,] / (median_reads[1,] - median_reads[7,] ) )*100), #hairpin miRNAs
                     ((median_reads[5,] + median_reads[4,])/ (median_reads[1,]  - median_reads[7,] ) )*100, # average of tRNAs
                     ((median_reads[6,] / (median_reads[1,]  - median_reads[7,]) )*100), #snoRNA reads
                     # ((median_reads[7,] / median_reads[1,])*100), #rRNA reads
                     ((median_reads[8,] / (median_reads[1,] - median_reads[7,]))*100) #other ncRNA reads
                     )
  
  perc_reads$cat = c("Mature miRNAs","Hairpin miRNAs","tRNA Reads","snoRNA Reads",
                     # "rRNA Reads",
                     "other ncRNAs") #rownames(perc_reads)
  View(perc_reads)
  
  data2plot = melt(perc_reads) #colnames(data2plot)
  colnames(data2plot) = c("RNAtype" , "timepoint", "value")
  
  
  if(select == 1){
    data2plot <- subset(data2plot, RNAtype %in% c("tRNA Reads"))
  }
  View(data2plot)
  
  datapiplot <- data2plot[which(data2plot$timepoint == "d0"),]
  
  de<-data.frame("Remaining","d0", 100 - sum(datapiplot$value))
  names(de) = names(datapiplot)
  datapiplot = rbind(datapiplot, de)
  
  # Add label position
  datapiplot <- datapiplot %>%
    arrange(desc(RNAtype)) %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  datapiplot
  
  #------------------------------------------------------------------
  q <- ggplot(datapiplot, aes(x = "", y = value, fill = RNAtype)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_label_repel(aes(y = lab.ypos,
                         label = paste0(format(round(value, 1), nsmall = 1), "%")),
                     color = "black", nudge_x = 0.5, nudge_y = 0.5
                     , show.legend = FALSE)+
    # scale_fill_manual(values = rainbow(7)) +
    # scale_fill_brewer(palette = "Pastel1") +
    scale_colour_brewer(hcl.colors(5,"ArmyRose"))+# colorBlindBlack8) +
    ggtitle(title) +
    theme_void(base_size = 14,
               base_family = "Helvetica") #+
    #theme(text=element_text(size=6,  family="Helvetica"))

  #------------------------------------------------------------------
  # Stack plot to explore each timepoint
  p <- ggplot(data2plot, aes(fill=RNAtype, y=value, x=timepoint)) + 
          geom_bar(position="stack", stat="identity", width=0.5) +
      scale_fill_viridis(discrete = T) +
      # geom_text(aes(label = format(round(value, 1), nsmall = 1)),
      #         position=position_stack(0.5), vjust=-0.5, size=2.5,
      #         color = "black",fontface = "bold", check_overlap = T) +
      ggtitle(title) +
      theme_ipsum() + xlab("Timepoints") + ylab("Percentage of total reads")
  
  return(q) # returns the pie chart
  # return(p)    #returns the stacked bar plot
}

# Percent stacked Plot 
# ggplot(data2, aes(fill=RNAtype, y=value, x=timepoint)) + 
#   geom_bar(position="fill", stat="identity") +
#   # scale_fill_viridis(discrete = T) +
#   geom_text(aes(label = format(round(value, 2), nsmall = 2)), 
#             position=position_fill(1.0), vjust=-0.5, size=2.5, 
#             color = "black",fontface = "bold", check_overlap = T) +
#   ggtitle("Hematopoetic Cells") + theme_ipsum() +  xlab("Timepoints")

save_png <-function(a, a_trna, filenamestring){
  
  save <- ggpubr::ggarrange(
    a, a_trna, 
    common.legend = FALSE, legend = "right",
    ncol = 2 , nrow = 1
  )
  
  #solving saving as pdf: https://github.com/thomasp85/ggraph/issues/152
  ggsave(
    filename = paste0(filenamestring, "_ncRNAtype_distributions_","1bcde",".4.png"),
    save, width=9.38, height=3.96  # width=8, height=7.4
  )
}


a <- snc_RNAdis_plot(file = "CM.annotation.report.csv", title = "Cardiomyocytes", cell_type="CM",select = 0);a
a_trna <- snc_RNAdis_plot(file = "CM.annotation.report.csv", title = "Cardiomyocytes", cell_type="CM",select = 1);a_trna
save_png(a, a_trna, filenamestring = "CM")

b <- snc_RNAdis_plot(file = "FB.annotation.report.csv", title = "Fibroblasts", cell_type="FB", select=0);b
b_trna <- snc_RNAdis_plot(file = "FB.annotation.report.csv", title = "Fibroblasts", cell_type="FB", select=1);b_trna
save_png(b, b_trna, filenamestring = "FB")


c <- snc_RNAdis_plot(file = "EC.annotation.report.csv", title = "Endothelial Cells", cell_type="EC", select=0);c
c_trna <- snc_RNAdis_plot(file = "EC.annotation.report.csv", title = "Endothelial Cells", cell_type="EC", select=1);c_trna
save_png(c, c_trna, filenamestring = "EC")

d <- snc_RNAdis_plot(file = "HC.annotation.report.csv", title = "Hematopoietic Cells", cell_type="HC", select=0);d
d_trna <- snc_RNAdis_plot(file = "HC.annotation.report.csv", title = "Hematopoietic Cells", cell_type="HC", select=1);d_trna
save_png(d, d_trna, filenamestring = "HC")

# saving the pie chart for the AMI paper
ggpubr::ggarrange(
  a, b, c, d, #align='hv',
  # labels = c("b", "c", "d", "e"),
  common.legend = TRUE, legend = "right",
  ncol = 4 , nrow = 1
)

ggsave(
  filename = paste0("ncRNAtype_distributions_","1bcde",".4.pdf"),
  # plot = last_plot(),  
  width=9.38, height=3.96  # width=8, height=7.4
)




