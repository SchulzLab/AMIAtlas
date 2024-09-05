library(ggplot2)
suppressWarnings(suppressMessages(library("DESeq2")))
library(data.table)
library(ggfortify)
library(dplyr)

exp_dir = "./expression/"
filename="_mature_normalized_CPM.1.txt"


cltyps = c("CM", "EC", "FB", "HC")

pr <- function(ct){
	mir_expr <- read.csv2(paste0(exp_dir,
	                    #  ct,
	                     "/", ct,
	                     filename),
	              header = TRUE, sep = "\t")

	print(head(mir_expr))
	print(colnames(mir_expr))

	rownames(mir_expr) <- mir_expr$miRNA
	mir_expr <- mir_expr[, -1]
	#print(head(mir_expr))
	mir_expr <- as.data.frame(sapply(mir_expr, as.numeric))

	t_mir_expr <- transpose(mir_expr)

	# get row and colnames in order
	colnames(t_mir_expr) <- rownames(mir_expr)
	rownames(t_mir_expr) <- colnames(mir_expr)
	t_mir_expr$timepoint=factor(rep(c("0","1","3","7","14","28"),each=4),
		levels = c("0","1","3","7","14","28")	)
	t_mir_expr$replicate=factor(rep(c("r1","r2", "r3","r4"),6))
	print(str(t_mir_expr))

	# df <- t_mir_expr[seq(1:24)]			#iris[c(1,2,3,4)]

	return(t_mir_expr)

}

# --------------------CM: plot PCA -----------------------
ct = "CM"
t_mir_expr <- pr(ct = "CM")
df <- t_mir_expr[seq(1:24)]

# plotpca <- function(){

	pdf(paste0("PCA_smallRNA_",ct,".pdf"), 7,7)

	theme_set(theme_minimal())
	autoplot(prcomp(df, scale. = TRUE), 
		data = t_mir_expr, 
		colour = 'replicate', 
		shape='timepoint', 
		size = 5) +
	  ggplot2::ggtitle(paste0(ct))

	dev.off()
# --------------------EC: plot PCA -----------------------
ct = "EC"
t_mir_expr <- pr(ct = "EC")
df <- t_mir_expr[seq(1:24)]

# plotpca <- function(){

	pdf(paste0("PCA_smallRNA_",ct,".pdf"), 7,7)

	theme_set(theme_minimal())
	autoplot(prcomp(df, scale. = TRUE), 
		data = t_mir_expr, 
		colour = 'replicate', 
		shape='timepoint', 
		size = 5) +
	  ggplot2::ggtitle(paste0(ct))

	dev.off()
# --------------------FB: plot PCA -----------------------

ct = "FB"
t_mir_expr <- pr(ct = "FB")
df <- t_mir_expr[seq(1:24)]

# plotpca <- function(){

	pdf(paste0("PCA_smallRNA_",ct,".pdf"), 7,7)

	theme_set(theme_minimal())
	autoplot(prcomp(df, scale. = TRUE), 
		data = t_mir_expr, 
		colour = 'replicate', 
		shape='timepoint', 
		size = 5) +
	  ggplot2::ggtitle(paste0(ct))

	dev.off()

# --------------------HC: plot PCA -----------------------

ct = "HC"
t_mir_expr <- pr(ct = "HC")
df <- t_mir_expr[seq(1:24)]

# plotpca <- function(){

	pdf(paste0("PCA_smallRNA_",ct,".pdf"), 7,7)

	theme_set(theme_minimal())
	autoplot(prcomp(df, scale. = TRUE), 
		data = t_mir_expr, 
		colour = 'replicate', 
		shape='timepoint', 
		size = 5) +
	  ggplot2::ggtitle(paste0(ct))

	dev.off()

