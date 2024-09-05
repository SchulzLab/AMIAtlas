# return the number of targets in cell type
load_library<-function(){
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    suppressMessages(library(ggplot2))
    library(ggrepel)
    library(tidyverse)
}
# this function calculates the number of targets for each miRNA
numtargets <- function(celltype){

    file = paste0(celltype,"_spearman_alltimepoints_anno.tsv")
    
    # check number of unique miRNAs
    cat("reading file ...", file)
    
    data<-read.csv2(file=paste0(dir,file)
                    , sep="\t"
                    ,header = TRUE) 
    
      data_filter <- data %>% 
                        dplyr::mutate(tmp.estimate = as.numeric(as.character(tmp.estimate))) %>%
                        filter(tmp.estimate < -0.4) %>%
                        count(mature)
    
    return(data_filter)    

}

dir = "./correlations/"
load_library()

dCM <- numtargets("CM")
dFB <- numtargets("FB")
dEC <- numtargets("EC")
dHC <- numtargets("HC")
nrow(dCM);nrow(dEC);nrow(dFB);nrow(dHC);


# merge dataframes from all cell types
dCM_EC <- merge(dCM, dEC, by="mature", all=TRUE); print(head(dCM))
dCM_EC_FB <- merge(dCM_EC, dFB, by="mature", all=TRUE)
colnames(dCM_EC_FB) <- c("mature", "nCM", "nEC", "nFB"); print(head(dCM_EC_FB))
dCM_EC_FB_HC <- merge(dCM_EC_FB, dHC, by="mature", all=TRUE)
colnames(dCM_EC_FB_HC) <- c("mature", "nCM", "nEC", "nFB", "nHC"); #print(head(dCM_EC_FB_HC))


# compute average per cell type
nrow(dCM_EC)
View(dCM_EC_FB_HC)

dCM_EC_FB_HC$restCM <- rowMeans(dCM_EC_FB_HC[ , c("nEC","nFB","nHC")], na.rm=TRUE)
dCM_EC_FB_HC$restEC <- rowMeans(dCM_EC_FB_HC[ , c("nCM","nFB","nHC")], na.rm=TRUE)
dCM_EC_FB_HC$restFB <- rowMeans(dCM_EC_FB_HC[ , c("nEC","nCM","nHC")], na.rm=TRUE)
dCM_EC_FB_HC$restHC <- rowMeans(dCM_EC_FB_HC[ , c("nCM","nEC","nFB")], na.rm=TRUE)

dCM_EC_FB_HC$sum <- rowSums(dCM_EC_FB_HC[ , c("nCM","nEC","nFB", "nHC")], na.rm=TRUE)

dCM_EC_FB_HC <- dCM_EC_FB_HC %>% mutate(rankCM = nCM/restCM
                                        ,rankEC = nEC/restEC
                                        ,rankFB = nFB/restFB
                                        ,rankHC = nHC/restHC) 

saveRDS(dCM_EC_FB_HC, file = "./results/dCM_EC_FB_HC.Rds")
dCM_EC_FB_HC <- readRDS(file = "./results/dCM_EC_FB_HC.Rds"); View(dCM_EC_FB_HC)



data <- dCM_EC_FB_HC %>% na.omit() #colnames(dCM_EC_FB_HC %>% na.omit() )

#%     ================== check expression ratio =============================== #

# get the miRNA expression files
test_expr_ratio_stats <-function(ct){
    
    exp_dir = "./expression/"
    filename="_mature_normalized_CPM.1.txt"
    
    cltyps = c("CM", "EC", "FB", "HC")
    ct = ct
    mir_expr <- read.csv2(paste0(exp_dir,
                     ct,
                     "/", ct,
                     filename),
              header = TRUE, sep = "\t")
    # View(mir_expr)
    mir_expr$ct = ct
    print( colnames(mir_expr) )
    colnames(mir_expr) = c("X", paste0(rep("D0_", 4), 1:4),
                            paste0(rep("D1_", 4), 1:4),
                            paste0(rep("D3_", 4), 1:4),
                            paste0(rep("D7_", 4), 1:4),
                            paste0(rep("D14_", 4), 1:4),
                            paste0(rep("D28_", 4), 1:4), "ct")
    
    # compute the ratio in the cell type of interest
    something <- mir_expr %>% 
                    pivot_longer( !c("ct", "X"),
                    names_to = "timepoint", 
                    values_to = "expr") %>% 
        rowwise() %>%
        dplyr::mutate(time = strsplit(timepoint,split='_', fixed=TRUE)[[1]][1]
                      ,logexpr = log2(as.numeric(expr) )) 
    
    # View(something)
    
    gene_expr_stats = something %>% 
        dplyr::group_by(time, X, ct) %>%
        dplyr::summarise(mean = mean(logexpr), sd = sd(logexpr))
    
    # View(gene_expr_stats)
    
    return(gene_expr_stats)

}
load_library()
a <- bind_rows(test_expr_ratio_stats(ct="CM"), test_expr_ratio_stats(ct="EC"))
a <- bind_rows(a, test_expr_ratio_stats(ct="FB"))
a <- bind_rows(a, test_expr_ratio_stats(ct="HC")) ; View(a)


# now find the max of each miRNA within each cell type
a_max = a %>% 
    dplyr::group_by(X, ct) %>%
    dplyr::summarise(mean = max(mean)); View(a_max)

# calculate the expression ratio
store=data.frame()

for (miR in (dCM_EC_FB_HC$mature)){
    
    # miR = "mmu-miR-465b-5p"
    a_max_sel <- a_max %>% filter(str_detect(X, paste0(miR, "\\b"))) 
    
    
    a_max_sel[1,4] = a_max_sel[1,3] / max(a_max_sel[2,3], a_max_sel[3,3], a_max_sel[4,3] )
    a_max_sel[2,4] = a_max_sel[2,3] / max(a_max_sel[1,3], a_max_sel[3,3], a_max_sel[4,3] )
    a_max_sel[3,4] = a_max_sel[3,3] / max(a_max_sel[1,3], a_max_sel[2,3], a_max_sel[4,3] )
    a_max_sel[4,4] = a_max_sel[4,3] / max(a_max_sel[1,3], a_max_sel[2,3], a_max_sel[3,3] )
    
    colnames(a_max_sel)[4] = "rel_exp_ratio"
    
    a_max_sel <- a_max_sel %>% pivot_wider(id_cols=!c("mean"),
                              names_from = "ct", 
                              values_from = "rel_exp_ratio") #%>% print()
    if(ncol(a_max_sel) > 4)
        store <- bind_rows(store, a_max_sel)
    
}

store <- merge(store, dCM_EC_FB_HC, by.x = "X", by.y = "mature", all=TRUE)

#================================================================================#

lollipop_plot_hc <- function(data){

    data = data %>% select(!c("NA")) %>%        # dCM_EC_FB_HC %>% 
                na.omit() %>%
                arrange(rankHC) %>% 
                filter((HC) > 1) %>%
                tail(5)
    data <- data %>% mutate(mature = paste0(gsub("mmu-", "", X)))
    data$mature <- factor(data$mature, levels = data$mature)
    
    # plot
    a <- ggplot(data, aes(x=mature, y=rankHC)) +
        geom_segment( aes(x=mature, xend=mature, y=8, yend=rankHC)) +
        geom_point( aes(size=(data$nHC/20)), color="black", fill=alpha("white", 0.0), alpha=0.7, shape=21, stroke=1) +
        scale_size(range = c(1,10), name="Size\n(unit)", breaks = c(1, 2.5, 5, 10)) +
        # geom_point( size=15, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1) +
        ylim(8, 17.5) +
        theme_minimal(base_size = 18) + 
        theme(text=element_text(color="black"),
              axis.text = element_text(color="black"))+
        ggtitle("Hematopoietic cells") +
        xlab("") + ylab("Ratio") +
        coord_flip()
    
    print(a)
    ggsave("./results/HC_specific_miRNAs.lollipop.2.pdf", height = 5.57, width = 5.57)
}
lollipop_plot_hc(store)

lollipop_plot_fb <- function(data){

    data = data %>%  select(!c("NA")) %>%             #dCM_EC_FB_HC %>% 
                na.omit() %>% 
                arrange(rankFB) %>% 
                filter(abs(FB) > 1, !(X %in% c("mmu-miR-6952-3p", "mmu-miR-7689-3p"))) %>%
                tail(5)
    data <- data %>% mutate(mature = paste0(gsub("mmu-", "", X)))
    data$mature <- factor(data$mature, levels = data$mature)
    
    # plot
    a <- ggplot(data, aes(x=mature, y=rankFB)) +
        geom_segment( aes(x=mature, xend=mature, y=3, yend=rankFB)) +
        geom_point( aes(size=(data$nFB/20)), color="black", fill=alpha("white", 0.0), alpha=0.7, shape=21, stroke=1) +
        scale_size(range = c(1,10), name="Size\n(unit)", breaks = c(1, 2.5, 5, 10)) +
        # geom_point( size=15, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1) +
        ylim(3, 15) +
        theme_minimal(base_size = 18) + 
        theme(text=element_text(color="black"),
              axis.text = element_text(color="black"))+
        ggtitle("Fibroblasts") +
        xlab("") + ylab("Ratio") +
        coord_flip()
    print(a)
    ggsave("./results/FB_specific_miRNAs.lollipop.1.pdf", height = 5.57, width = 5.57)
}
lollipop_plot_fb(store)

lollipop_plot_ec <- function(data){

    data = data %>% select(!c("NA")) %>%                           # dCM_EC_FB_HC %>% 
                na.omit() %>% 
                arrange(rankEC) %>% 
                filter(abs(EC) > 1) %>%
                tail(5)
    data <- data %>% mutate(mature = paste0(gsub("mmu-", "", X)))
    data$mature <- factor(data$mature, levels = data$mature)
    
    # plot
    a <- ggplot(data, aes(x=mature, y=rankEC)) +
        geom_segment( aes(x=mature, xend=mature, y=2.5, yend=rankEC)) +
        geom_point( aes(size=(data$nEC/20)), color="black", fill=alpha("white", 0.0), alpha=0.7, shape=21, stroke=1) + #fill=alpha("#FAD5A5", 0.3),
        scale_size(range = c(1,10), name="Size\n(unit)", breaks = c(1, 2.5, 5, 10)) +
        # geom_point( size=15, color="red",  alpha=0.7, shape=21, stroke=1) +
        ylim(2.5, 10) +
        theme_minimal(base_size = 18) + 
        theme(text=element_text(color="black"),
              axis.text = element_text(color="black"))+
        ggtitle("Endothelial cells") +
        xlab("") + ylab("Ratio") +
        coord_flip()
    print(a)
    ggsave("./results/EC_specific_miRNAs.lollipop.1.pdf", height = 5.57, width = 5.57)
}
lollipop_plot_ec(store)

lollipop_plot_cm <- function(data){

    data = data %>% select(!c("NA")) %>%                                               #dCM_EC_FB_HC %>% 
                na.omit() %>%
                arrange(rankCM) %>% 
                filter(abs(CM) > 1) %>%
            tail(5)
    data <- data %>% mutate(mature = paste0(gsub("mmu-", "", X)))
    data$mature <- factor(data$mature, levels = data$mature)
    
    # plot
    a <- ggplot(data, aes(x=mature, y=rankCM)) +
        geom_segment( aes(x=mature, xend=mature, y=2.5, yend=rankCM)) +
        
        geom_point( aes(size=(data$nCM/20)), color="black", fill=alpha("white", 0.0), alpha=0.7, shape=21, stroke=1) +
        scale_size(range = c(1,10), name="Size\n(unit)", breaks = c(1, 2.5, 5, 10)) +
        # geom_point( size=15, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1) +
        ylim(2.5, 10) +
        theme_minimal(base_size = 18) + 
        theme(text=element_text(color="black"),
              axis.text = element_text(color="black"))+
        ggtitle("Cardiomyocytes") +
        xlab("") + ylab("Ratio") +
        coord_flip()
        
    print(a)
    ggsave("./results/CM_specific_miRNAs.lollipop.1.pdf", height = 5.57, width = 5.57)
}
lollipop_plot_cm(store)

























