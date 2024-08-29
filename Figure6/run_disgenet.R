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



dotplot(xx2, showCategory = 20)
dotplot(ck, showCategory = 20)





ck1 <- compareCluster(geneCluster = something
                     , fun = enrichDO
                     # , OrgDb = 'org.Hs.eg.db'
)
dotplot(ck1)



x = enrichDO(list,                                 # DOSE 
             ont           = "DO",
             pvalueCutoff  = 0.05,
             pAdjustMethod = "BH",
             minGSSize     = 5,
             maxGSSize     = 500,
             qvalueCutoff  = 0.05,
             readable      = FALSE);                                       
bubble_plot(x, pattern="", filelabel=paste0(celltype, "_DOSE_allfunctions"), height = 20, width = 10)
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
pdf(paste0(celltype,"_DOSE_treeplot.pdf"), width = 15)
treeplot(edox2)
dev.off()

dgn <- enrichDGN(list);                               # DisGeNET enrichment
bubble_plot(dgn, pattern="", filelabel=paste0(celltype, "_Disegnet_allfunctions"), height = 120, width = 20)
edox <- setReadable(dgn, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
pdf(paste0("/Volumes/Elements_1/smallRNA_MI/disease_assoc/", celltype,"_treeplot.pdf"), width = 15)
treeplot(edox2)
dev.off()

ncg <- enrichNCG(list)                                # NetworkCancer
bubble_plot(ncg, pattern="", filelabel=paste0(celltype,"_NetworkCancer_allfunctions"), height = 120, width = 10)

# ================================== =================================== #
# Can we draw a plot to show the biological concepts enriched 
# among the ceRNA genes (mouse)  

loadLibraries()
# get the gene symbols of ceRNA network
celist <- c("Actb", "Adamts5", "Adamtsl3", "Akap12", "Bicc1"
,"Camk1d", "Cd84", "Col12a1", "Col1a1","Col4a3","Col5a3","Col8a1","Crispld2"
,"Ddr2","Dio2","Enc1","Egfr","Fbn1","Fn1","Gfpt2","Igf1"
,"Ipcef1","Iqgap2","Kirrel","Lamc1","Lox","Mmp19","Myo5a"
,"Pax5","Pcsk6","Pdgfra","Plekho2","Plxnc1","Prex1"
,"Rassf4","Reps2","Robo1","Runx1"
,"Serpine2","Sgms2","Sh3pxd2b","Slc16a3","Spi1","Sulf1"
,"Tcf21","Tifab")
# convert to ENTREZ ID
up_names = gprofiler2::gconvert(celist
                    , organism = "mmusculus" #"hsapiens",
                    ,target = "ENTREZGENE_ACC" # "ENSG"
); #View(up_names)

edo_bp <- enrichGO(up_names$target        #de
                , ont="BP"
                , 'org.Mm.eg.db')   # enrichDGN(de)
edo_bp <- pairwise_termsim(edo_bp)
edox <- setReadable(edo_bp, 'org.Mm.eg.db', 'ENTREZID')
categorys <- c("extracellular matrix organization" 	
               ,"bone development"
               , "positive regulation of fibroblast proliferation"
               # ,"urogenital system development"
               ,"extracellular matrix disassembly"
               )
a <- mutate(edo_bp, qscore = -log(p.adjust, base=10))
barplot(a, x="qscore", showCategory = categorys, width = 1.5, title = "Biological process", font.size = 14)    
ggsave("GOBP_enrichfunctions_barplots_plot.pdf", width=8, height = 2.8)  


# a1 <- dotplot(edo_bp, showCategory=categorys) + ggtitle("Biological process")
p1 <- cnetplot(edox
               , showCategory = categorys #6 #categorys
               , categorySize="p.adjust"
               ,'gem'
               , foldChange=up_names$name)+ 
    ggtitle("Biological process"); View(edox@result)



edo_mf <- enrichGO(up_names$target        #de
                , ont="MF"
                , 'org.Mm.eg.db')   # enrichDGN(de)
edo_mf <- pairwise_termsim(edo_mf)
categorys <- c("extracellular matrix structural constituent"
               ,"glycosaminoglycan binding"
               ,"growth factor binding"
               ,"integrin binding"
)
b <- mutate(edo_mf, qscore = -log(p.adjust, base=10)) 
barplot(b, x="qscore",  width = 1.5, showCategory = categorys, title = "Molecular function", font.size = 14) 
ggsave("GOMF_enrichfunctions_barplots_plot.pdf", width=8, height = 2.8)   


edox <- setReadable(edo_mf, 'org.Mm.eg.db', 'ENTREZID') ## convert gene ID to Symbol
p2 <- cnetplot(edox
               , showCategory = categorys    #5
               , categorySize="p.adjust"
               ,'gem'
               , foldChange=up_names$name) + 
    ggtitle("Molecular function")
p2;View(edox@result)


edo_cc <- enrichGO(up_names$target        #de
                   , ont="CC"
                   , 'org.Mm.eg.db')   # enrichDGN(de)
edo_cc1 <- pairwise_termsim(edo_cc)
edox <- setReadable(edo_cc1, 'org.Mm.eg.db', 'ENTREZID')   ## convert gene ID to Symbol
categorys <- c("collagen-containing extracellular matrix"
               # ,"collagen trimer"
               # ,"basement membrane"
               ,"fibrillar collagen trimer"
               # ,"cell leading edge"
               # ,"membrane raft"
                )
c <- mutate(edo_cc1, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", width = 0.3, showCategory = categorys) 
ggsave("GOCC_enrichfunctions_barplots_plot.pdf"
       # , width=8
       # , height = 1.8
       )
p3 <- cnetplot(edox
               , showCategory = categorys    #7
               , categorySize="p.adjust"
               ,'gem'
               ,foldChange=up_names$name) + 
    ggtitle("Cellular component")
View(edox@result)

p <- cowplot::plot_grid(a, b, c
                        ,nrow=3          
                        ,rel_heights = c(.8, .8, .7)
                        ,rel_widths = c(.8, .8, .8));  #labels=LETTERS[1:3]
ggsave("GO_enrichfunctions_barplots_plot.pdf", width=8, height = 6)

p <- cowplot::plot_grid(p1, p2, p3, ncol=3
                        # , labels=LETTERS[1:3]
                        , rel_widths=c(.8, .8, 0.8))
ggsave("cnet_plot.pdf", width=20, height = 6)

# ==================================================================#

# library
library(ggplot2)
library(viridis)
library(hrbrthemes)

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
    theme_ipsum() + 
    geom_bar(position="dodge", stat="identity")

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Graph
ggplot(data, aes(fill=condition, y=value, x=condition)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("Studying 4 species..") +
    facet_wrap(~specie) +
    theme_ipsum() +
    theme(legend.position="none") +
    xlab("")



# heatplot(edox, showCategory=5)

#specific selections for better view
bubble_plot(x, pattern="lymph|imm", filelabel="lymph")
bubble_plot(dgn, pattern="hepa|immune", filelabel="hepa")



