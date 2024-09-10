
## ========================= load libraries ===================== 
loadlibraries <- function(){
  library("ggrepel")
  library("ggplot2")
  library("viridis")
  library("hrbrthemes")
  library("matrixStats")
  library("reshape2")
  library("readr")
  library("dplyr")    # Data manipulation
  library("ggpubr")
  library("readxl")
}


#To Do: Check for MDS plot for choosing the samples, getting rid of the outlier samples
# Outlier samples: ctrl_HC_3, ctrl_FB_2, d28_EC_3

snc_RNAdis_plot <- function(title, cell_type, select){
  

  annotation_report <- read_xlsx(paste0(
                                    ".",
                                    "small_RNA-report.xlsx"), 
                                    sheet=paste0(cell_type, "_smallRNA"), 
                                    range="A1:O25")

  #choose the columns
  # [1]"Sample name(s)" [2]"Total Input Reads" [3]"Trimmed Reads (all)"  [4]"Trimmed Reads (unique)" [5]"All miRNA Reads" [6]"Filtered miRNA Reads" [7]"Unique miRNAs"         
  # [8]"Hairpin miRNAs" [9]"mature tRNA Reads" [10]"primary tRNA Reads" [11]"snoRNA Reads" [12]"rRNA Reads" [13]"ncRNA others" [14]"mRNA Reads" [15] "Remaining Reads"
  
  EC <- annotation_report[, c(1, 2, 6, 8, 9, 10, 11, 12, 13, 15)] #head(EC)

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

  
  data2plot = melt(perc_reads) #colnames(data2plot)
  colnames(data2plot) = c("RNAtype" , "timepoint", "value")  
  
  if(select == 1){
    data2plot <- subset(data2plot, RNAtype %in% c("tRNA Reads"))
  }
  
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

a <- snc_RNAdis_plot(title = "Cardiomyocytes", cell_type="CM",select = 0);a

b <- snc_RNAdis_plot(title = "Fibroblasts", cell_type="FB", select=0);b

c <- snc_RNAdis_plot(title = "Endothelial Cells", cell_type="EC", select=0);c

d <- snc_RNAdis_plot(title = "Hematopoietic Cells", cell_type="HC", select=0);d

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
  width=9.38, height=3.96
)
