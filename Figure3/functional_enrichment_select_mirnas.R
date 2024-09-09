# this is functional enrichment of the targets of the miRNAs in each cluster
# =========================== libraries =========================
library(maSigPro)
library(pheatmap)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gprofiler2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ReactomePA))
# =========================== set the options ====================

library("optparse")
 
option_list = list(
    make_option(c("-c", "--ct"), type="character", default=NULL, 
              help="cell type, CM, EC, FB or HC", metavar="character"),

    make_option(c("-k", "--klusters"), type="character", default=NULL, 
              help="number of clusters you are interested in", metavar="character")     #,

    # make_option(c("-r", "--corr"), type="character", default=NULL, 
    #         help="correlation cutoff, default: -0.4", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$ct)){
  print_help(opt_parser)
  stop("provide cell type, CM, EC, FB or HC", call.=FALSE)
} else if (is.null(opt$k)){
  print_help(opt_parser)
  stop("provide number of clusters you are interested in", call.=FALSE)
}
#  else if (is.null(opt$corr)){
#   print_help(opt_parser)
#   stop("provide correlation cutoff", call.=FALSE)
# }

# =============================== functions ============================ #
get_targtes <- function(
        corr_table,          #dataframe with corr coeff and p values
        corr_cutoff,           # correlation cutoff
        p_val_cutoff,          #p-value cutoff
        mir){
    
    corr_table_sel = corr_table %>% 
        filter(tmp.p.value < p_val_cutoff, tmp.estimate < corr_cutoff , mature == mir)     

    return(corr_table_sel$gene)
}

# generate the table of functions
func_enrichment_targets_clusterProfiler_table <- function(
        directory, 
        celltype, 
        corr_cutoff, 
        mirnas,
        cluster){
    
    print("here in function <func_enrichment_targets_clusterProfiler_table>")
    
    #read corr coeff and select targets 
    print(paste0(directory, celltype,"_spearman_alltimepoints_anno.tsv"))
    corr_df = read.csv(paste0(directory, celltype,"_spearman_alltimepoints_anno.tsv"), sep="\t") # output FROM spearman.R #glimpse(corr_df)
    print(head(corr_df, 3))     # check to see if read the correct file
    
    i=1
    merge_all =  data.frame()
    for (mir_arg in mirnas) { 
        
        
        cat("Iteration: ",i,"/",length(mirnas),": ",mir_arg,"\n")
        sel_target_genes_3p = get_targtes(corr_table= corr_df, 
                                          corr_cutoff = corr_cutoff, 
                                          p_val_cutoff = 0.05, 
                                          mir = c(paste0(mir_arg, ""))) 
        
        eg3 = bitr(str_to_title(sel_target_genes_3p), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db"); 
        
        bp2 = data.frame()
        # ======================= Compute GO functional enrichment ===========================
        ego3 <- enrichGO(gene = eg3$ENTREZID,
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'ENTREZID',
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         readable = TRUE,
                         pvalueCutoff  = 0.01 #,
                         # qvalueCutoff  = 0.05
        )
        
        if(!is.null(ego3)) {
            ck <- setReadable(ego3, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
        }else{
            ck = NULL
        }
        print("---------------------- Enriching target genes with GO Functions ------------------------")
              
        if(!is.null(ck)){
            
            if(nrow(ck) > 0){
                
                bp <- pairwise_termsim(ck)
                # bp2 <- simplify(bp, cutoff=0.3, by="p.adjust", select_fun=min)
                  
                
                # bp3 <- simplify(bp, cutoff=0.8, by="p.adjust", select_fun=min)
                
                
                bp2 <- simplify(bp, cutoff=0.1, by="p.adjust", select_fun=min)  
                
                
            }
        }
        
        kk3 <- enrichKEGG(gene = eg3$ENTREZID,
                          organism     = 'mmu',
                          pvalueCutoff = 0.05)
        if(!is.null(kk3)){
            ck <- setReadable(kk3, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")
        }else{
            ck <- NULL
        }
        print("---------------------- Enriching target genes with KEGG PATHWAYs  ------------------------")
        
        if(!is.null(ck)){
            if(nrow(ck)>0){
                
                kk3 = as.data.frame(ck)
                kk3$ONTOLOGY = "KEGG"
                bp2 <- bind_rows(as.data.frame(bp2), kk3) 
            }
        }
        
        print("---------------------- Enriching target genes with WIKI Pathways ------------------------")
        # x <- enrichWP(gene = eg3$ENTREZID, organism = "Mus musculus") 
        # if(nrow(as.data.frame(x)) > 0) {     #!is.null(x)){
        #     xr <- as.data.frame(setReadable(x, OrgDb = "org.Mm.eg.db", keyType="ENTREZID"))
        #     xr$ONTOLOGY = "WP"
        #     # print(head(xr))   
        #     # write.table(xr, paste0("tabs/", celltype, "_",mir_arg,"_clusterProfiler_WP_func_enrich_results_3p.2.tsv"), sep = "\t", quote=FALSE, row.names=FALSE)
        #     bp2 <- bind_rows(as.data.frame(bp2), as.data.frame(xr))
        # }        
        
        print("---------------------- Enriching target genes with REACTOME PATHWAYs ------------------------")
        x <- enrichPathway(gene=eg3$ENTREZID, pvalueCutoff = 0.05, readable=TRUE, organism = "mouse")
        if(nrow(as.data.frame(x)) > 0){ #!is.null(x)){
            #       
            # print(head(x,2))
            x = as.data.frame(x)
            x$ONTOLOGY = "RP"
            # write.table(x, paste0("tabs/", celltype, "_",mir_arg,"_clusterProfiler_RP_func_enrich_results_3p.2.tsv"), sep = "\t", quote=FALSE, row.names=FALSE)
            bp2 <- bind_rows(as.data.frame(bp2), x)
            
        }        
        
        i = i+1
        merge_all = bind_rows(merge_all, as.data.frame(bp2))
    }
    
    #create directory structure to store results
    dir.create("results/")

    merge_all = merge_all[order(merge_all$p.adjust, decreasing = FALSE), ]  # order by adjusted p.value
    merge_all$CLUSTER = paste0("m",opt$ct,cluster)                          # write the cluster name
    
    write.table(merge_all, 
                paste0(
                    "./results/", 
                    celltype, "_Cluster",cluster,"_func_enrich_",corr_cutoff,".tsv"), 
                    append = FALSE, sep = "\t", quote=FALSE, row.names=FALSE)
}


# -----------------  step 1: get the clustering results again  ----------------- #
dir = "/Volumes/Elements_1/smallRNA_MI/miRBase_counts/"  # directory to read the sample information from
celltype = "FB"
k = 6 #specify number of clusters desired

if(celltype == "CM"){
    samples_CM <- read.table(file.path(paste0(dir, "CM/samples_CM_new.txt")),
                         header = TRUE)
}else if(celltype == "EC"){
    samples_CM <- read.table(file.path(paste0(dir, "EC/samples_EC_remout.txt")), header=TRUE);

}else if(celltype == "HC"){
    samples_CM <- read.table(file.path(paste0(dir,"HC/samples_HC_remout.txt")), header=TRUE)
}else if(celltype == "FB"){
  samples_CM <- read.table(file.path(paste0(dir,"FB/samples_FB_remout.txt")), header=TRUE)
}else{
    stop("input celltype parameter not correct!!")
}


if(celltype == "CM"){
    rownames(samples_CM) = c(
        paste("ctrl", rep("CM", 4), 1:4, sep = "_"),
        paste("d1", rep("CM", 4), 1:4, sep = "_"),
        paste("d3", rep("CM", 4), 1:4, sep = "_"),
        paste("d7", rep("CM", 4), 1:4, sep = "_"),
        paste("d14", rep("CM", 4), 1:4, sep = "_"),
        paste("d28", rep("CM", 4), 1:4, sep = "_")
    )     #head(samples_CM,5)
}else if(celltype == "EC"){
    rownames(samples_CM) = c(paste("ctrl", rep("EC",4), 1:4, sep="_"), 
         paste("d1", rep("EC",4), 1:4, sep="_"),
         paste("d3", rep("EC",4), 1:4, sep="_"),
         paste("d7", rep("EC",4), 1:4, sep="_"),
         paste("d14", rep("EC",4), 1:4, sep="_"),
         paste("d28", rep("EC",2), 1:2, sep="_"),
         paste("d28", rep("EC",1), 4, sep="_"))
}else if(celltype == "HC"){
    rownames(samples_CM) = c(paste("ctrl", rep("HC",2), 1:2, sep="_"), 
         paste("ctrl", rep("HC",1), 4, sep="_"),
         paste("d1", rep("HC",4), 1:4, sep="_"),
         paste("d3", rep("HC",4), 1:4, sep="_"),
         paste("d7", rep("HC",4), 1:4, sep="_"),
         paste("d14", rep("HC",4), 1:4, sep="_"),
         paste("d28", rep("HC",4), 1:4, sep="_"))
}else if(celltype == "FB"){
  rownames(samples_CM) =  c(paste("ctrl", rep("FB",1), 1, sep="_"), 
                            paste("ctrl", rep("FB",2), 3:4, sep="_"),
                            paste("d1", rep("FB",4), 1:4, sep="_"),
                            paste("d3", rep("FB",4), 1:4, sep="_"),
                            paste("d7", rep("FB",4), 1:4, sep="_"),
                            paste("d14", rep("FB",4), 1:4, sep="_"),
                            paste("d28", rep("FB",4), 1:4, sep="_"))
}

if(celltype == "CM"){
    data.mir.CM = read.table(file.path(paste0(dir, "CM/CM_mature_normalized_CPM.txt")), 
                             header = TRUE); 
    head(data.mir.CM, 5)
}else if(celltype == "EC"){
    data.mir.CM = read.table(file.path(paste0(dir, "EC/EC_mature_normalized_CPM.txt")), 
                             header=TRUE)
    print(head(data.mir.CM, 5))
    drops <- c("d28_EC_3.mature")
    data.mir.CM = data.mir.CM[ , !(colnames(data.mir.CM) %in% drops)]
    print(head(data.mir.CM,3))
}else if(celltype == "HC"){
    data.mir.CM = read.table(file.path(paste0(dir,"HC/HC_mature_normalized_CPM.txt")), header=TRUE)
    drops <- c("ctrl_HC_3")
    data.mir.CM = data.mir.CM[ , !(colnames(data.mir.CM) %in% drops)]
    print(head(data.mir.CM,3))
    
}else if(celltype == "FB"){
  data.mir.CM = read.table(file.path(paste0(dir,"FB/FB_mature_normalized_CPM.txt")), header=TRUE)
  drops <- c("ctrl_FB_2.mature")
  data.mir.CM = data.mir.CM[ , !(colnames(data.mir.CM) %in% drops)]
  print(head(data.mir.CM,3))
  
}

design <- make.design.matrix(samples_CM, degree = 5)
if((celltype != "HC") & (celltype != "FB")){
    rownames(data.mir.CM) <- data.mir.CM[, 1]
    data.mir.CM[, 1] <- NULL
}

colnames(data.mir.CM) = rownames(samples_CM)
data.mir.CMn = data.mir.CM[rowSums(data.mir.CM) > 10, ]

cat("Dimensions before setting threshold: ", dim(data.mir.CM), "\n")
cat("Dimensions after setting threshold: ", dim(data.mir.CMn), "\n")

# running maSigpro for significance analysis
fitCM <-
    p.vector(
        data = data.mir.CMn,
        design,
        Q = 0.05,
        MT.adjust = "BH",
        min.obs = 20
    )

# recompute miRNA clusters
res <- pheatmap(
    fitCM$SELEC,
    scale = "row",
    clustering_distance_rows = "correlation",
    cluster_cols = F,
    show_rownames = F,
    border_color = NA,
    cutree_rows = k,
    filename = paste0("Heatmap_significant_miR_",celltype,"_mod.2.pdf")
)

# get the miRNAs in each cluster
CM.clust <- data.frame(cbind(fitCM$SELEC,
                             cluster = cutree(res$tree_row,
                                              k = k)))
for (x in 1:k) {
    rowname_vec <- rownames(CM.clust[CM.clust$cluster == x, ])
    cat("Running Cluster ", x, ": ", length(rowname_vec), "\n")
    
    if(x>5){
      func_enrichment_targets_clusterProfiler_table(
        directory = "./correlations/",
        celltype = celltype,
        corr_cutoff = -0.4,
        mirnas = rowname_vec,   # mirnas$x
        cluster  = x
      )
    }
    
}
