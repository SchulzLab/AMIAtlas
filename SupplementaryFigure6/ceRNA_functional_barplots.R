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