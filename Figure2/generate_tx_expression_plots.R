# This program takes a gene input and plots (in log2 scale) the mean expression along the timepoints
# for all the celltypes
# argument 1 = ensembl gene id
# argument 2 = gene name (as appears in the title of the plot)
# argument 3 = name appended to the file name

#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
	make_option(c("-e", "--ens"), type="character", default=NULL, 
              help="ensembl gene id", metavar="character"),

	make_option(c("-g", "--gname"), type="character", default=NULL, 
              help="gene name (as appears in the title of the plot)", metavar="character"),

	make_option(c("-o", "--outname"), type="character", default="out.txt", 
              help="name appended to the output file name", metavar="character")

	make_option(c("-c", "--celltype"), type="character", default=NULL, 
        help="cell type, one of CM, EC, FB, HC", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$ens)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as the ensemble gene id", call.=FALSE)
} else if (is.null(opt$g)){
  print_help(opt_parser)
  stop("provide gene name (as appears in the title of the plot)", call.=FALSE)
} else if (is.null(opt$o)){
	stop("name appended to the output file name, preferable the gene name / symbol", call.=FALSE)
}else if (is.null(opt$c)){
	stop("provide interested cell type", call.=FALSE)
}



suppressPackageStartupMessages(library(dplyr, tidyverse))
library(tidyr)
library(ggplot2)
library(egg)
library(hrbrthemes)
library(stringr)
library(gridExtra)
library(cowplot)



plot_tx_expression <- function(mtr_stats_wider_all, gene_name){

	p <- ggplot(mtr_stats_wider_all, aes(x=factor(time) )) #color = str_wrap(rowname, 10)		
	
	p <- p + geom_line( aes(y=mean , group = ct, color = ct), linewidth = 1.4) +
				geom_point(position= position_dodge2(0.1), aes(y=mean, color = ct), size=1, alpha = 0.8) +
				geom_errorbar(aes(ymin=mean - sd, color = ct, ymax=mean + sd),
	                  					position=position_dodge(0.05), width=.2 , alpha = 0.8) + 
				scale_color_manual(values=c("#CC6666", "#9999CC"))

	p <- p + theme_bw() +
  	# scale_y_continuous(breaks = scales::breaks_pretty(n = 0)) +
    theme(legend.position="bottom",
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
          legend.margin=margin(t=-10),
          # legend.key.size =  unit(0.3, "in"),
          legend.text=element_text(size=5),
          text = element_text(size=18, color = "black"),
		  axis.text.x = element_text(color = "black"),
		  axis.text.y = element_text(color = "black"),
          axis.title = element_text(face="bold"),
          axis.title.x.bottom = element_text(face="bold", size=8),
          axis.title.y = element_text(face="bold", size=8),
          axis.title.y.right = element_text(face="bold", size=8)) +
    labs(	x = "Timepoints", colour = NULL,
    		y = "Mean Normalized Expression (Log2 TPM)", 
    		title = paste0(gene_name, "")) +   # : gene expression")) +
    guides(color=guide_legend(nrow=1,byrow=TRUE))

    print(p)
	
}

dir = "./expression/"
cltyps = opt$c				#c("CM") # "CM", "EC", "FB", "HC")
gene_expr_all = data.frame()
for (ct in head(cltyps,4)) {

	cat("Reading for the celltype: ", ct, "\n")
	cat("reading from file: ", paste0(dir, "norm_gene_counts_",ct,".csv"))
	gene_expr = read.csv(paste0(dir, "norm_gene_counts_",ct,".csv"),
					header = TRUE, sep = ",") %>% filter(str_detect(X, paste0(opt$ens, '.'))) 
	

	#modify the colnames
	gene_expr$ct = ct

	if(ct == "FB"){
		# missing "d3_FB_2"  "d3_FB_3"
		gene_expr$d3_FB_2 = gene_expr$d3_FB_1
		gene_expr$d3_FB_3 = gene_expr$d3_FB_4
		print(colnames(gene_expr))

		# rearrange columns
		col_order <- c("X", paste0(rep("d0_FB_", 4), 1:4),
									paste0(rep("d1_FB_", 4), 1:4),
									paste0(rep("d3_FB_", 4), 1:4),
									paste0(rep("d7_FB_", 4), 1:4),
									paste0(rep("d14_FB_", 4), 1:4),
									paste0(rep("d28_FB_", 4), 1:4), "ct")
		gene_expr <- gene_expr[, col_order]
		print("reordered----------------")
		print(head(gene_expr))


	}

	if(ct == "HC"){

		print(colnames(gene_expr))
		gene_expr <- gene_expr %>% mutate(d28_HC_4 = mean(c(d28_HC_1, d28_HC_2, d28_HC_3)))
		print(head(gene_expr))
		
		# # rearrange columns
		col_order <- c("X", paste0(rep("d0_HC_", 4), 1:4),
									paste0(rep("d1_HC_", 4), 1:4),
									paste0(rep("d3_HC_", 4), 1:4),
									paste0(rep("d7_HC_", 4), 1:4),
									paste0(rep("d14_HC_", 4), 1:4),
									paste0(rep("d28_HC_", 4), 1:4), "ct")
		gene_expr <- gene_expr[, col_order]
	
		print("reordered----------------")
		print(head(gene_expr))

	}

	colnames(gene_expr) = c("X", paste0(rep("D0_", 4), 1:4),
								paste0(rep("D1_", 4), 1:4),
								paste0(rep("D3_", 4), 1:4),
								paste0(rep("D7_", 4), 1:4),
								paste0(rep("D14_", 4), 1:4),
								paste0(rep("D28_", 4), 1:4), "ct")
	print(head(gene_expr))

	# append to dataframe to take in other celltypes
	gene_expr_all <- bind_rows(gene_expr_all, gene_expr)
	print(head(gene_expr_all))

}

something <- gene_expr_all %>% 
	    				pivot_longer(!c("X", "ct"), names_to = "timepoint", values_to = "expr") %>% 
	    				rowwise() %>%
	    				dplyr::mutate(time = strsplit(timepoint,split='_', fixed=TRUE)[[1]][1]
	                  					,logexpr = log2(as.numeric(expr) )) 


  	gene_expr_stats = something %>% 
						dplyr::group_by(time, X, ct) %>%
						dplyr::summarise(mean = mean(logexpr), sd = sd(logexpr))

	print(gene_expr_stats$time)
	gene_expr_stats$time <- factor(gene_expr_stats$time, levels=c("D0", "D1", "D3","D7","D14" ,"D28" ))

	print(gene_expr_stats)

	plot_tx_expression(gene_expr_stats, gene_name=opt$gname) #"Stk39")

	ggsave(paste0("./results/",opt$outname,"_gene_expression_log2.1.pdf"), 
			units = "in", width = 7, height = 3)

#"Stk39
