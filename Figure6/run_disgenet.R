loadLibraries <- function(){
    suppressMessages(library(magrittr) )
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    suppressMessages(library(DOSE))
    suppressMessages(library(ggplot2))
    suppressMessages(library(forcats))
    suppressMessages(library(enrichplot))
    suppressMessages(library(data.table))
    suppressMessages(library(stringr))
    suppressMessages(library(clusterProfiler))
    library(gprofiler2)
}


library(AnnotationDbi)
library(org.Mm.eg.db)
data(geneList, package="DOSE")

str(geneList)
AnnotationDbi::select(org.Mm.eg.db
                      # , keys=genes
                      , columns='ENTREZID', keytype='SYMBOL')

all <- read.csv2("https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt", sep="\t", header = TRUE) %>% as.data.frame()
View(all)

#example dataset
data(geneList); #View(geneList)
de = names(geneList)[1:100]

bubble_plot <- function(x, pattern, filelabel, height, width){
    
    y <- as.data.frame(x) %>%
        mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
    
    if(pattern == ""){
        y_select = y
    }else{
        y_select <- y %>% filter(str_detect(Description, regex(pattern, ignore_case = TRUE))) 
    }
    
    View(y); # for preview of the entire list #View(y_select)
    ggplot(y_select, showCategory = 20, 
           aes(richFactor, fct_reorder(Description, richFactor))) + 
            geom_segment(aes(xend=0, yend = Description)) +
            geom_point(aes(color=p.adjust, size = Count)) +
            scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
            scale_size_continuous(range=c(2, 10)) +
            theme_minimal() + 
            xlab("rich factor") +
            ylab(NULL) + 
            ggtitle("Enriched Disease Ontology")
    
    ggsave(filename = paste0("/Volumes/Elements_1/smallRNA_MI/disease_assoc/try_",filelabel,".pdf"), 
           width = width,
           height = height, limitsize = FALSE)
    
}


# this function returns the list of ENTREZ gene ids (max differentiating at a timepoint),
# given a celltype. This function thus returns the list of the genes for each cell type
return_orth_genelist <- function(celltype){
    
    # read then maSigPro summary file; concatenate and unique
    summ_genes <- read.csv2(paste0("/Volumes/Elements_1/AG_Simon/CATS_trimming/masigpro/", celltype,"_summary_genes.txt"), sep=",") ; #View(summ_genes)
    
    # find human orthologous genes
    if(celltype == "HC"){
        list = summ_genes$Group.1vsGroup.0
    }else{
        list = summ_genes$Group.3vsGroup.0
    }
    
    human_genes = gorth(list,                           # summ_genes$Group.1vsGroup.0, 
                        target_organism = "hsapiens", 
                        source_organism = "mmusculus",
                        numeric_ns = "ENTREZGENE_ACC"); 
    
    #find the orthogonal genes of all mouse genes
    all_mouse_genes_tohuman = gorth(as.integer(all$X5..genome.build) %>% unique(), # summ_genes$Group.1vsGroup.0, 
                        target_organism = "hsapiens", 
                        source_organism = "mmusculus",
                        numeric_ns = "ENTREZGENE_ACC");
    
    universe <- gconvert(all_mouse_genes_tohuman$ortholog_name,
                         organism = "hsapiens",
                         target = "ENTREZGENE_ACC")
    
    # find their entrez gene IDs
    up_names = gconvert(human_genes$ortholog_name,
                        organism = "hsapiens",
                        target = "ENTREZGENE_ACC" # "ENSG"
    ); #View(up_names)
    
    # make a condensed list
    something = list()
    something$interest = up_names$target
    something$unniverse = universe$target
    return(something)
}

loadLibraries()
something = list()
str(return_orth_genelist(celltype = "CM"))
something$CM <- unlist(return_orth_genelist(celltype = "CM")$interest)
something$EC <- unlist(return_orth_genelist(celltype = "EC")$interest)
something$FB <- unlist(return_orth_genelist(celltype = "FB")$interest)
something$HC <- unlist(return_orth_genelist(celltype = "HC")$interest)
universse = unlist(return_orth_genelist(celltype = "CM")$unniverse)

str(something)

# enrichDGN1 = enrichDGN(gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, 
#                        minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, readable = FALSE)

enrichDGN1 <- function(gene, universe = universse){
    return(enrichDGN(gene, universe = universe))
}

library(org.Hs.eg.db)
ck <- compareCluster(geneCluster = something
                     , fun = 'enrichDGN1'   
                     # , OrgDb = 'org.Hs.eg.db'
)
# ck <- setReadable(ck, OrgDb = 'org.Hs.eg.db', keyType="ENTREZID")
xx2 <- pairwise_termsim(ck)

categorys <- c("Abnormal coordination", 	
               "Abnormal endocardium morphology"
               , "acute aortic dissection"
               ,"Aortic Aneurysm"
               ,"Congenital muscular dystrophy (disorder)"
               ,"Restrictive cardiomyopathy"
              , "Myocardial Failure"
              ,"Atrial Septal Defects"
              ,"Myocardial infarction, stroke"
              ,"Myocardial Reperfusion Injury"
              ,"Ischemic cardiomyopathy"
              ,"Hypertrophic obstructive cardiomyopathy"
              ,"Acute myocardial ischemia"
              ,"Aortic Valve Stenosis"
              ,"Dilatation of aorta"
              ,"Muscular Dystrophy"
              ,"Hypoxia",
              "Respiratory Failure"
               )
View(ck@compareClusterResult)
pdf("/Volumes/Elements_1/smallRNA_MI/disease_assoc/compared_disgenet.3.pdf", height = 4, width = 8)
dotplot(xx2, showCategory = categorys, font.size=14) #10) #
dev.off()









