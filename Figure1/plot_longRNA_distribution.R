loadLibraries <-function() {
    
    suppressMessages(library(ggplot2))
    suppressMessages(library(viridis))
    # suppressMessages(library(hrbrthemes))
    library(ggrepel)
    library(readxl)
    library(ggplot2)
    library(tidyverse)
    cat("loaded libraries\n")   
    
}


plot_longRNA_disrtibution <- function(longRNA_FB, title, ct) {
    pc <- longRNA_FB %>% filter(biotype == "PC") %>% select(median) %>% as.list()
    val_pc <- median(pc$median)
    
    ri <- longRNA_FB %>% filter(biotype == "retained introns") %>% select(median) %>% as.list()
    val_ri <- median(ri$median)
    
    lnc <- longRNA_FB %>% filter(biotype == "lncRNAs") %>% select(median) %>% as.list()
    val_lnc <- median(lnc$median)
    
    pse <- longRNA_FB %>% filter(biotype == "pseudogenes") %>% select(median) %>% as.list()
    val_pse <- median(pse$median)
    
    pt <- longRNA_FB %>% filter(biotype == "processed transcripts") %>% select(median) %>% as.list()
    val_pt <- median(pt$median)
    
    s <- sum(median(pc$median), median(ri$median), median(lnc$median), median(pse$median), median(pt$median))
    
    datapiplot <- data.frame(
        class = c("PCGs", "Retained introns", "lncRNAs", "Pseudogenes", "processed transcripts"),
        n = c(median(pc$median), median(ri$median), median(lnc$median), median(pse$median), median(pt$median)),
        prop = c((val_pc/s)*100,(val_ri/s)*100, (val_lnc/s)*100, (val_pse/s)*100, (val_pt/s)*100)
    )
    
    print(datapiplot)
    colorBlindBlack8  <- c(#"#000000", 
                           "#E69F00", "#56B4E9", "#009E73", 
                           "#F0E442", "#0072B2"
                           #, "#D55E00", "#CC79A7"
                           )
    datapiplot <- datapiplot %>%
        arrange(desc(class)) %>%
        mutate(lab.ypos = cumsum(prop) - 0.5*prop)
    
    p <- ggplot(datapiplot, aes(x = "", y = prop, fill = class)) +
        geom_bar(width = 1, stat = "identity", color = "white") +
        coord_polar("y", start = 0)+
        geom_text(aes(y = lab.ypos, label=""
                      # label = format(round(prop, 1), nsmall = 1)
        ), 
        color = "black", 
        nudge_x = 0.0, 
        nudge_y = 5) +
        scale_colour_brewer(hcl.colors(5,"ArmyRose"))+      # colorBlindBlack8) +
        ggtitle(title) +
        theme_void(base_size = 14,
                   base_family = "Helvetica")
    
    print(p)        # see the plot in rstudio editor
    return(p)    
    
}

hcl.colors(5,"Green-Orange")
loadLibraries()
longRNA_FB1<- read_xlsx("./ncRNA_quantification_logs.xlsx", 
                       sheet = "FB", range = "A1:C37")
a <- plot_longRNA_disrtibution(longRNA_FB1, title="Fibroblast", ct="FB")
longRNA_EC<- read_xlsx("./ncRNA_quantification_logs.xlsx", 
                        sheet = "EC", range = "A1:C37")
b <- plot_longRNA_disrtibution(longRNA_EC, title="Endothelial cells", ct="EC")

longRNA_CM<- read_xlsx("./ncRNA_quantification_logs.xlsx", 
                       sheet = "CM_new", range = "A1:C37")
c <- plot_longRNA_disrtibution(longRNA_CM, title="Cardiomyodytes", ct="CM")

longRNA_HC <- read_xlsx("./ncRNA_quantification_logs.xlsx", 
                       sheet = "HC_new", range = "A1:C37")
d <- plot_longRNA_disrtibution(longRNA_HC, title="Hematopoietic cells", ct="HC")



ggpubr::ggarrange(
    a, b, c, d, #align='hv', # labels = c("b", "c", "d", "e"),    
    common.legend = TRUE, legend = "right",
    ncol = 4 , nrow = 1
)

ggsave(
    filename = paste0(paste0("./",
                             "AMI_longRNA_distribution.1.pdf")),
    # plot = last_plot(),  
    width=9.38, height=3.96
)







