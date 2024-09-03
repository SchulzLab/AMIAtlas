## ============================= on HC and EC long RNA seq data ==============================
setwd("/Users/ranjan/time_course_analysis/masigpro/")
library(readr)
library(maSigPro)
library(MASS)

library("optparse")

option_list = list(
    make_option(c("-c", "--celltype"), type="character", default=NULL, 
              help="cell type of interest, either CM, EC, FB or HC", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$c)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as the cell type", call.=FALSE)
} 
# ======================================= install.packages("UpSetR") =============================================
library("UpSetR")
pdf(paste0("upset_",opt$c,".2.pdf"),
        width = 5, 
        height = 6,
        onefile=F)
listInput <- list(d1vsd0 = get1$summary$Group1vsGroup0, 
                  d3vsd0 = get1$summary$Group2vsGroup0, 
                  d7vsd0 = get1$summary$Group3vsGroup0,
                  d14vsd0 = get1$summary$Group4vsGroup0,
                  d28vsd0 = get1$summary$Group5vsGroup0)
upset(fromList(listInput), 
                order.by = "freq",
                text.scale = c(1.5, 1.5, 1.5, 0.7, 2, 0))
#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
dev.off()



