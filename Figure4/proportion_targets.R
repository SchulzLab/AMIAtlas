# this is the script for figure 2 where we find the overlap of miRNA target 
# genes, from each miRNA cluster,  with the 
# STEM clusting output of various gene cluster profiles, for each cell type
# loading all functions
loadLibraries <- function() {
  library(maSigPro)
  suppressMessages(library(reshape))
  library(ggplot2)
  library(pheatmap)
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(gprofiler2))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(AnnotationDbi))
  suppressMessages(library(org.Mm.eg.db))
  suppressMessages(library(gprofiler2))
}


#1. get the miRNA clusters from the heatmap
get_mir_clusters <- function(celltype, k, corr_coeff) {

  cat("entering function get_mir_clusters\n")

  dir <- "./expression/"
  if (celltype == "CM") {
    samples <- read.table(file.path(paste0(dir, "CM/samples_CM_new.txt")),
                          header = TRUE)
  }else {
    samples <- read.table(file.path(paste0(dir,
                                          #  celltype,
                                           "/samples_",
                                           celltype,
                                           "_remout.txt")),
                          header = TRUE)
  }
  print(samples)
  print(rownames(samples))
	
  if (celltype == "CM")
    rownames(samples) <- c(paste("ctrl", rep(celltype,4), 1:4, sep = "_"),
                           paste("d1", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d3", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d7", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d14", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d28", rep(celltype, 4), 1:4, sep = "_"))
  if (celltype == "EC")
    rownames(samples) <- c(paste("ctrl", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d1", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d3", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d7", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d14", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d28", rep(celltype,3), 1:3, sep="_"))

  if (celltype == "FB" | celltype == "HC")
    rownames(samples) <- c(paste("ctrl", rep(celltype, 3), 1:3, sep = "_"),
                           paste("d1", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d3", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d7", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d14", rep(celltype, 4), 1:4, sep = "_"),
                           paste("d28", rep(celltype, 4), 1:4, sep = "_"))

  data.mir.CM <- read.table(file.path(paste0(dir,
                                            #  celltype,
                                             "/", celltype,
                                             "_mature_normalized_CPM.1.txt")),
                            header = TRUE)
  print(head(data.mir.CM, 5))
	
  design <- make.design.matrix(samples, degree = 5)
  rownames(data.mir.CM) <- data.mir.CM[,1]
  data.mir.CM[, 1] <- NULL
  colnames(data.mir.CM) <- rownames(samples)
  data.mir.CMn <- data.mir.CM[rowSums(data.mir.CM) > 10, ]
     
  fitCM <- p.vector(data = data.mir.CMn, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
	
	print(fitCM$i)
	
  #get the miRNAs in each cluster
  res <- pheatmap(fitCM$SELEC,
                  scale = "row",
                  clustering_distance_rows = "correlation",
			            cluster_cols = F,
                  show_rownames = F,
                  border_color = NA,
                  cutree_rows = k,
                  filename = paste0("./results/Heatmap_significant_miR_", celltype, ".", 
                                    k,
                                    "_mod.2.pdf")
                  )

	CM.clust <- data.frame(cbind(fitCM$SELEC, 
								cluster = cutree(res$tree_row, 
												k = k)))

  #put the targets of the miRNAs in the cluster as a list
	df <- list()
  for (x in 1:k) {
    rowname_vec <- rownames(CM.clust[CM.clust$cluster==x,])
    cat("Cluster ", x, ": ",length(rowname_vec),"\n")
    print(rowname_vec)
    
    temp <- get_mircluster_targets(rowname_vec, celltype, corr_coeff) %>%
            dplyr::select(gene) %>%
            as.list() #%>% unlist()
    
    df[[x]] <- temp$gene
  # write.table(df, paste0("cluster", x, celltype,"_targets.csv"))
    cat("\n")
  }

  return(df)

}


# get the targets of the miRNAs in each of the clusters, provided as list 
get_mircluster_targets <- function(mir_list, celltype, corr_coeff){

  cat("get gene targets of ", mir_list, "\n")

  dir = "./correlations/"

	#read the correlations
  file = paste0(dir, celltype, "_spearman_alltimepoints_anno.tsv")
    
  cat("reading file: ", file,"\n")
  data <- read.csv2(file=paste0(file)
                    , sep = "\t"
                    , header = TRUE) 

  data_filter <- data %>% 
                      dplyr::mutate(tmp.estimate = as.numeric(as.character(tmp.estimate))) %>%
                      dplyr::filter(tmp.estimate < corr_coeff) %>%
                      dplyr::filter(mature %in% mir_list) #== "mmu-let-7a-2-3p")

  return(data_filter)	

}

# find overlap, calculate ratio and compute fishers exact test
find_overlap <- function(df, gene_profiles, celltype){

  cat("in function find_overlap")
  #print(str(gene_profiles));
  cat("--------------------------------")
  #print((df[[1]])) #print((gene_profiles[[1]]))
  cat("--------------------------------\n")
  #print((df[[2]])) #print((gene_profiles[[2]]))
  #initialize a dataframe
  a <- data.frame(x = rep(paste(substring(celltype, 1, 1), 1:length(df))
                        , length(gene_profiles))
                  , y = rep(paste("X", 1:length(gene_profiles)), length(df))
                  , value = runif(1)
                  , pval = runif(1)
                  , intersect = runif(1)
                  , logpval = runif(1)
                  , stringsAsFactors = FALSE)
  count <- 1

  totalmiRtargets = c()
	for (j in 1:length(df)){	
			
			totalmiRtargets = append(totalmiRtargets, unique(df[[j]]))
	}

	totalmiRtargets = unique(totalmiRtargets)

	totalgenes = c()
	for (j in 1:length(gene_profiles)){	
			
			# totalgenes = totalgenes + length(unique(gene_profiles[[j]]))
			totalgenes = append(totalgenes, unique(gene_profiles[[j]]))

	}

  for (i in 1:length(gene_profiles)){

    for (j in 1:length(df)){

      #if from masigpro, convert to ensem gene symbols
					tar <- gconvert(gene_profiles[[i]]
						, organism="mmusculus"
						, target = "ENSG", numeric_ns = "ENTREZGENE_ACC")
					
					num <- length(unique(intersect(toupper(unique(df[[j]]))
											,toupper(tar$name	)				#gene_profiles[[i]])
											))
								)

					den = length(tar$name)							#gene_profiles[[i]])

					a$x[count] <- paste0(substring(celltype, 1, 1),j) #; print(a$x[count] )
					a$y[count] <- paste0("X",i)

					print(num)
					print(den)
					print(num/den)
					numb = (num/den)
					denb = length(totalmiRtargets)/length(totalgenes)

					a$value[count] = numb #/ denb #(length(totalmiRtargets)/length(totalgenes))

					S = length (intersect(toupper(setdiff(totalmiRtargets, unique(df[[j]])))
										,toupper(tar$name)))			#gene_profiles[[i]])))

					tar1 <- gconvert(setdiff(totalgenes, gene_profiles[[i]]), organism="mmusculus"
						, target = "ENSG", numeric_ns = "ENTREZGENE_ACC")$name


					G <- length(intersect(toupper(tar1)   #tar$name)
										, toupper(setdiff(totalmiRtargets , unique(df[[j]])))))


					#compute fisher's exact
					dat <- data.frame(
						"targets" = c(num, length(df[[j]])-num),
						"not-targets" = c(S, G), 
						row.names = c("mirT", "Non-mirT"),
						stringsAsFactors = FALSE
					)
					colnames(dat) <- c("ingeneprofile", "notingeneprofile")

					# print(dat) #print(test); print(test$p.value)
					test <- fisher.test(dat)			
					
					if (round(test$p.value,3) == 1){
						a$pval[count] = 0 #test$p.value
						a$logpval[count] = 0 
					}
					else{
						a$pval[count] = test$p.value
						a$logpval[count] = -log10(test$p.value)
					}

					# if(-log2(test$p.value) < 4.32)
					# 	a$logpval[count] = 0 
					
					a$intersect[count] = num
					count = count+1
			}
	}

	pvalues <- a$pval	
	a$padj <- p.adjust(pvalues, method="BH") 			#"bonferroni")
	a$logpadj <- -log2(a$padj)

	# > p.adjust(pvalues,method="hochberg")
	# > p.adjust(pvalues,method="BH")

	a <- a %>% mutate(modlogpadj = case_when(logpadj < 4.32 ~ 0
											, (logpadj > 4.32 & logpadj != Inf) ~ logpadj
											, (logpadj == Inf) ~ 0))
	
	print(a)
	return(a)	

}

# function to plot the ratio and the fishers exact test p value
plot_heatmap <- function(d, celltype, width, height, corr_coeff, k){

  d <- d %>% mutate(y = forcats::fct_relevel(y, paste0("X", 1:10)))

  ggplot(d, aes(y
          ,forcats::fct_rev(x)
          ,fill = value
          , size = modlogpadj #logpadj #logpval
        )) +
    geom_point(shape = 21, stroke = 0) +
    geom_hline(yintercept = seq(.5, 10.5, 1), size = .2) +
    geom_vline(xintercept = seq(.5, 15.5, 1), size = .2) +
    scale_x_discrete(position = "top") +
    scale_size(range = c(0, 12)) +
    scale_fill_gradient(low = "orange",
                        high = "blue",
                      breaks = c(0, 0.1, (max(d$value)))
                      , labels = c("0", "", as.character(round(max(d$value),2)))    #"0.2")
                      , limits = c(0,  max(d$value))) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", #"top", 
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                               label.position = "top",
                               title.position = "top",#"right", 
                               order = 1),
           fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
    labs(size = "-log(adjp.val)", fill = "Normalized \nratio", x = NULL, y = NULL)
	
	ggsave(paste0("./results/Norm_intersect.masigpro.",celltype,".", k, "." ,corr_coeff,".2.pdf")
			, width = width, height = height)  # )

}

read_masigpro_results <- function(celltype, k){
	print("inside function to read_masigpro_results ...")
	print(getwd())

	if(celltype=="CM" | celltype=="EC")
		d <- read.csv2(paste0("./significant_genes/",celltype,"_significant_genes_anno.txt")
					, sep=",")
	else
		d <- read.csv2(paste0("./significant_genes/",celltype,"_significant_genes_remout_anno.txt")
					, sep=",")
	
  if (celltype == "CM")
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
							# filename=paste0("gene_clusters/",celltype,"_gene_cluster_heatmap.pdf"), 
							# width = 3.5, 
							# height =  8							
							) 

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

	# Scale the data 
	# scaled_data <- t(scale(t(d_select)))    #scale_rows(d_select)		#scale(mtcars_data) 

	return(something)
	
}


loadLibraries()
#------------------------------------------------------------------#
# this section generates the dot plots 
#------------------------------------------------------------------#
for (corr_coeff in seq(-0.4,-0.9, -0.1)){
	df <- get_mir_clusters(celltype="CM", k=5, corr_coeff)
	# read the gene clusters from STEM output
	# gene_profiles <- read_STEM_profiles(celltype="CM")
	#generate the gene clusters from HCL on masigpro results
	gene_profiles <- read_masigpro_results(celltype="CM", k=10)
	dataf <- find_overlap(df, gene_profiles, celltype="CM")
	plot_heatmap(dataf, celltype="CM", width = 6.24, height = 5.30, corr_coeff, k=5)	#8

}

for (corr_coeff in seq(-0.4,-0.9, -0.1)){
	df <- get_mir_clusters(celltype="EC", k=5, corr_coeff)
	# gene_profiles <- read_STEM_profiles(celltype="EC")
	gene_profiles <- read_masigpro_results(celltype="EC", k=10)
	dataf <- find_overlap(df, gene_profiles, celltype="EC")
	plot_heatmap(dataf, celltype="EC", width = 8.24, height = 5.30, corr_coeff, k=5)			#10
}

for (corr_coeff in seq(-0.4,-0.9, -0.1)){
	df <- get_mir_clusters(celltype="FB", k=6, corr_coeff)
	# gene_profiles <- read_STEM_profiles(celltype="FB")
	gene_profiles <- read_masigpro_results(celltype="FB", k=10)
	dataf <- find_overlap(df, gene_profiles, celltype="FB")
	plot_heatmap(dataf, celltype="FB", width = 8.24, height = 5.30, corr_coeff, k=6)			#12
}

for (corr_coeff in seq(-0.4,-0.9, -0.1)){
	df <- get_mir_clusters(celltype="HC", k=6, corr_coeff)
	# gene_profiles <- read_STEM_profiles(celltype="HC")
	gene_profiles <- read_masigpro_results(celltype="HC", k=10)
	dataf <- find_overlap(df, gene_profiles, celltype="HC")
	plot_heatmap(dataf, celltype="HC", width = 8.24, height = 5.30, corr_coeff, k=6)
}



#------------------------------------------------------------------#
# to compare the functions berween functions from the gene clusters and the 
# target genes from the miRNA clusters
#------------------------------------------------------------------#


#do functional enrichment of the list of genes and return the GO ids/terms
get_functions <- function(df,noconvert) {

  cat("here in the function to do functional enrichment")

  store <- list()

  for (j in 1:length(df)){
    
    #convert to gene id to entrez id
    if(noconvert==1){
      list <- gconvert(df[[j]],
                      organism = "mmusculus",
                      target = "ENTREZGENE_ACC",
                      numeric_ns = "",
                      mthreshold = Inf,
                      filter_na = TRUE
                    ) 
      gene_list = list$target
    } else{
      gene_list <- df[[j]]
    }    
   
    ggo <- enrichGO(gene = gene_list,
               OrgDb    = org.Mm.eg.db,
               ont      = "BP",
               pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
               readable = TRUE)

    store[[j]] <- ggo$ID    
  }
  return(store)
}

#this function finds the overlap of the functional annotations
find_overlap_goids <- function(df, gene_profiles, celltype){
  cat("in function find_overlap_goids:\n")

  print(length(df))
  print(length(gene_profiles))
  print(df[[1]])
  print(gene_profiles[[1]])

  #datastructure to store the results
  a <- data.frame(x = rep(paste(substring(celltype, 1, 1), 1:length(df))
                        , length(gene_profiles))
                  , y = rep(paste("X", 1:length(gene_profiles)), length(df))
                  , value = runif(1)
                  , pval = runif(1)
                  , intersect = runif(1)
                  , logpval = runif(1)
                  , stringsAsFactors = FALSE)

  count <- 1
  #total functions covered in the miRtarget list
  totalmiRtargets = c()
	for (j in 1:length(df)){			
    totalmiRtargets = append(totalmiRtargets, unique(df[[j]]))
	}
	totalmiRtargets = unique(totalmiRtargets)

  #total functions covered in the gene cluster list
	totalgenes = c()
	for (j in 1:length(gene_profiles)){			
    totalgenes = append(totalgenes, unique(gene_profiles[[j]]))
	}
  print("here...")
  for (i in 1:length(gene_profiles)){

    for (j in 1:length(df)){
					
      num <- length(unique(intersect(toupper(unique(df[[j]]))
                                    ,toupper(gene_profiles[[i]]	)				#gene_profiles[[i]])
                                    )
                          )
                    )
      print(num)
     

					den = length(gene_profiles[[i]])							#gene_profiles[[i]])

					a$x[count] <- paste0(substring(celltype, 1, 1),j) #; print(a$x[count] )
					a$y[count] <- paste0("X",i)

					print(num)
					print(den)
					print(num/den)
					numb = (num/den)
					denb = length(totalmiRtargets)/length(totalgenes)

					a$value[count] = numb #/ denb #(length(totalmiRtargets)/length(totalgenes))
					
					S = length (intersect(toupper(setdiff(totalmiRtargets, unique(df[[j]])))
										,toupper(gene_profiles[[i]])))			#gene_profiles[[i]])))

					tar1 <- setdiff(totalgenes, gene_profiles[[i]])

					G <- length(intersect(toupper(tar1)   #tar$name)
										, toupper(setdiff(totalmiRtargets , unique(df[[j]])))))


					#compute fisher's exact
					dat <- data.frame(
						"targets" = c(num, length(df[[j]])-num),
						"not-targets" = c(S, G), 
						row.names = c("mirT", "Non-mirT"),
						stringsAsFactors = FALSE
					)
					colnames(dat) <- c("ingeneprofile", "notingeneprofile")

					# print(dat) #print(test); print(test$p.value)
					test <- fisher.test(dat)			
					
					if (round(test$p.value,3) == 1){
						a$pval[count] = 0 #test$p.value
						a$logpval[count] = 0 
					}
					else{
						a$pval[count] = test$p.value
						a$logpval[count] = -log10(test$p.value)
					}

					# if(-log2(test$p.value) < 4.32)
					# 	a$logpval[count] = 0 
					
					a$intersect[count] = num
					count = count+1
			}
	}

	pvalues <- a$pval	
	a$padj <- p.adjust(pvalues, method="BH") 			#"bonferroni")
	a$logpadj <- -log2(a$padj)

	# > p.adjust(pvalues,method="hochberg")
	# > p.adjust(pvalues,method="BH")

	a <- a %>% mutate(modlogpadj = case_when(logpadj < 4.32 ~ 0
											, (logpadj > 4.32 & logpadj != Inf) ~ logpadj
											, (logpadj == Inf) ~ 0))
	
	return(a)	

}




#-----------------------------------------------------------------------------#
# this function computes the functional enrichment for each miRNA, and returns table of functions
# table of functions
# Note: this function works on one miRNA at a time
#-----------------------------------------------------------------------------#
functions_permiRNA <- function(miR, celltype, coef){
  cat("here in function functions_permiRNA\n")  
  a <- get_mircluster_targets(miR, celltype, coef)  

  # write.table(a, 
  #             file=paste0("temp.",celltype, ".", miR,".functions.txt"),
  #             quote = FALSE,
  #             row.name = FALSE)
  cat("before conversion:", length(a$gene))
  #convert ot entrex gene ids
  list <- gconvert(a$gene,
                      organism = "mmusculus",
                      target = "ENTREZGENE_ACC",
                      numeric_ns = "",
                      mthreshold = Inf,
                      filter_na = TRUE
                    ) 
  
  gene_list = list$target
  cat("computing functional enrichment for ",
            (length(gene_list)),"targets\n")
  ggo <- enrichGO(gene = gene_list,
              OrgDb    = org.Mm.eg.db,
              ont      = "ALL", #c("MF", "BP"),
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.05,
              readable = TRUE)

  ggo %>% arrange(desc(p.adjust)) %>% head() %>% print() 
  ggo <- ggo %>% arrange((p.adjust))
  write.table(ggo, 
              file=paste0(celltype, ".", miR,".functions.txt"),
              sep = "\t",
              quote = FALSE,
              row.name = FALSE)
}