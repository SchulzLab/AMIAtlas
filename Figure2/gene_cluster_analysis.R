# this is the script for fig 2 to do the 

loadLibraries <- function(){

    suppressMessages(library(reshape))			# install.packages("reshape")
    suppressMessages(library(ggplot2))
    suppressMessages(library(pheatmap))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyr))
    suppressMessages(library(enrichplot))
    suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
    suppressMessages(library("tximport"))
    suppressMessages(library("RColorBrewer"))
    suppressMessages(library("grid"))
    suppressMessages(library("ggplotify"))
    suppressMessages(library("DESeq2"))

    # if(!require("xlsx")) {
    #     install.packages("xlsx")
    #     library("xlsx")
        
    # } else{
    #     library("xlsx")
    # }

    if (!require("maSigPro")) {
        BiocManager::install("maSigPro")
        suppressMessages(library(maSigPro))
    } else
        suppressMessages(library(maSigPro))
    if (!require("UpSetR")){
        install.packages("UpSetR")
        suppressMessages(library("UpSetR"))
    }else
        suppressMessages(library("UpSetR"))
    suppressMessages(library(clusterProfiler))
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    print("loaded all libraries ...")
}

#====================================================================================
# setwd("/Volumes/Elements/AG_Simon/CATS_trimming/")
# print("hello: Running differential gene expression analysis for CM samples")
# paste("Current working directory:", getwd())
#====================================================================================

read_masigpro_results <- function(celltype, k){
	print("inside function to read_masigpro_results ...")
	print(getwd())

	if(celltype=="CM" | celltype=="EC")
		d <- read.csv2(paste0("gene_clusters/",celltype,"_significant_genes_anno.txt")
					, sep=",")
	else
		d <- read.csv2(paste0("gene_clusters/",celltype,"_significant_genes_remout_anno.txt")
					, sep=",")
	
	if(celltype=="CM")
		d_select <- d %>% 
					dplyr::select(-Gene.name, -Gene.Synonym, -s) %>% 					# -id) 
					distinct()															# keep id column and remove duplicate expressions
	
	if(celltype=="EC" | celltype=="FB" | celltype=="HC")
		d_select <- d %>% 
					dplyr::select(-Gene.name, -Gene.Synonym) %>% 					# -id) 
					distinct()

	rownames(d_select) = d_select$id
		
	d_select <- d_select %>% dplyr::select(-id) 
	d_select <- mutate_all(d_select, function(x) as.numeric(as.character(x)))
	
    # res <- pheatmap(d_select,
    #     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    #         scale="row",
    #         clustering_distance_rows = "correlation",
    #         cluster_cols=F,
    #         show_rownames=F,
    #         border_color=NA,
    #         treeheight_row = 0,
    #         cutree_rows = 10,
    #         filename=paste0("gene_clusters/",celltype,"_Heatmap_significant_mod1.1.pdf"),
    #         width = 3.5,
    #         height  = 8
    # )

	clustered_data <- pheatmap(d_select, 								
							cluster_rows = TRUE, 
							cluster_cols=FALSE,
							show_colnames = TRUE, 
							main = "Customized Heatmap", 
							 scale="row",
							clustering_distance_rows = "correlation",
							color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), #c("#FFEDA0", "#F03B20"))(100), 
							fontsize_row = 8, 
							fontsize_col = 10, 
							cellwidth = 25, 
							cellheight = 15, 
							border_color = NA, 
							kmeans_k = k, 
							legend = TRUE, 
							legend_breaks = seq(-4, 4, by = 1.0), 
							legend_labels = c("-4", "-3", "-2", "-1",  "0",  "1",  "2" , "3",  "4" ),   
							fontsize_legend = 12,
							filename=paste0("gene_clusters/",celltype,"_gene_cluster_heatmap.pdf"), 
							width = 3.5, 
							height =  8							
							) 

	
	# d_select.clust <- data.frame(cbind(d_select, 
	# 							cluster = cutree(res$tree_row, 
	# 											k = k)))

	# print(res$tree_row$order)  
	clustering_res <- as.data.frame(clustered_data$kmeans$cluster) 
	colnames(clustering_res) = "cluster_id"
	d_select.clust <- merge(d_select, clustering_res, by=0, all=TRUE)
	rownames(d_select.clust) = d_select.clust$Row.names

	# print(rownames(d_select.clust)) #print(head(d_select.clust))
	# put to list
	something = list()
	length(something) = k
	names(something) = paste0(rep("X"),1:k,"")

	for (x in 1:k) {
		rowname_vec <- rownames(d_select.clust[d_select.clust$cluster_id==x,])
		cat("Cluster ", x, ": ",length(rowname_vec),"\n")

		something[[x]] = rowname_vec
	}

	# scaled_data <- t(scale(t(d_select)))    #scale_rows(d_select)		#scale(mtcars_data) # Scale the data 
	return(something)
	
}

functional_analysis <- function(res, celltype){

	cat("in function functional_analysis ...")
	print(str(res))

	ck <- compareCluster(geneCluster = res
						, OrgDb = 'org.Mm.eg.db'
						, fun = enrichGO)
	ck <- setReadable(ck 
						, OrgDb = 'org.Mm.eg.db'
						, keyType="ENTREZID")
	
	write.table(ck
				, paste0("gene_clusters/",celltype,"_functional_enrichments.before.txt"), sep="\t")

	xx2 <- pairwise_termsim(ck, method="JC")

	y <- as.data.frame(xx2) %>% rowwise() %>%
					mutate(GeneRatio =  Count / as.numeric(sub("\\d+/", "", GeneRatio))
							, logpadj =  -log10(p.adjust)
		) %>%
		filter(p.adjust < 0.05) #%>% select(Count, richFactor)
	
	
	# TODO: add column xlabel
	y %>% arrange(-GeneRatio, -logpadj) %>%
				write.table(paste0("gene_clusters/",celltype,"_functional_enrichments.1.txt"), 
				sep="\t", 
				row.names=FALSE)
    
	ggplot(data = y, aes(x = Cluster
						, y = Description
						, color = logpadj 				#`p.adjust`
						, size = GeneRatio)) + 
		geom_point() +
		scale_color_gradient(low = "red", high = "blue") +
		theme_bw() + 
		ylab("") + 
		xlab("") + 
		ggtitle(paste0("GO enrichment analysis of clusters in ",celltype))

	ggsave(paste0("gene_clusters/",celltype,"_gene_clusters_dotplot.pdf")
		,width=20
		,height=50
		, limitsize = FALSE)

	# TODO: select descriptions to show
	if(celltype=="CM")
		desc_list = c(
			"actin binding"
            # ,"vitamin binding"
            # ,"glycosaminoglycan binding"
            ,"GTPase regulator activity"
            # ,"kinase regulator activity"
            # ,"tubulin binding"
            , "extracellular matrix structural constituent"
            ,"structural constituent of ribosome"
            # ,"active transmembrane transporter activity"
            # ,"nucleoside-triphosphatase regulator activity"
            # ,"passive transmembrane transporter activity"
			)

	if(celltype=="EC")
		desc_list = c(
					# "structural constituent of ribosome"
					"extracellular matrix structural constituent"
					,"actin binding"
					,"ubiquitin-like protein transferase activity"
					# ,"active transmembrane transporter activity"
					,"GTPase regulator activity"
					, "cytoskeletal motor activity"
					,"nucleoside-triphosphatase regulator activity"
					# ,"glycosaminoglycan binding"
					# ,"GTPase regulator activity"
					# ,"passive transmembrane transporter activity"
				)
	if(celltype=="FB"){
		desc_list = c("structural constituent of ribosome"
						# ,"catalytic activity, acting on RNA"
						# ,"actin binding"
						,"GTPase regulator activity"
						# ,"nucleoside-triphosphatase regulator activity"
						# ,"passive transmembrane transporter activity"
						# ,"tubulin binding"
						# ,"protein serine/threonine kinase activity"
						, "heat shock protein binding"
						,"cell adhesion molecule binding"
						, "cell adhesion molecule binding"
						# ,"protein serine/threonine kinase activity"
		)}
	if(celltype=="HC"){
		desc_list = c(
			"actin binding"
			, "extracellular matrix structural constituent"
	# 		,"protein serine/threonine kinase activity"
			,"sulfur compound binding"
	# 		,"tubulin binding"
			,"GTP binding"
	# 		,"helicase activity"
			, "cytokine receptor binding"
			,"cytokine avtivity"
	# 		,"G protein-coupled receptor binding"
	# 		,"ubiquitin protein ligase activity",
	# 		"ribonucleoprotein complex binding",
	# 		"sulfur compound binding"
		)}
	
	y <- y %>% filter(Description %in% desc_list)

	ggplot(data = y, aes(x = Cluster
						, y = Description
						, color = logpadj 				#`p.adjust`
						, size = GeneRatio)) + 
		geom_point() +
		scale_color_gradient(low = "red", high = "blue") +
		theme_bw(base_size = 8) +
		ylab("") +  #ylab("GO functions") + 
		xlab("") + # xlab("Gene Clusters") + 
		theme(axis.text=element_text(size=18, colour = "black"))
		# ggtitle("GO enrichment analysis of clusters in CM (selected)")

    ggsave(paste0("gene_clusters/",celltype,"_gene_clusters_dotplot.selected.1.pdf"),
                    width = 12,
                    height = 2,
                    limitsize = TRUE)

}

loadLibraries()
set.seed(2024)
res <- read_masigpro_results(celltype="CM", k=10)
functional_analysis(res, celltype="CM")

res <- read_masigpro_results(celltype="EC", k=10)
functional_analysis(res, celltype="EC")

res <- read_masigpro_results(celltype="FB", k=10)
functional_analysis(res, celltype="FB")

res <- read_masigpro_results(celltype="HC", k=10)
functional_analysis(res, celltype="HC")



