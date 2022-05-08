#!/usr/bin/R

library(heatmaply)
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-r", "--rowmetadata"), type="character", default=NULL,
              help="metadata for rows", metavar="character"),
  make_option(c("-c", "--colmetadata"), type="character", default=NULL,
              help="metadata for cols", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (out file).\n", call.=FALSE)
}

mydata <- read.table(opt$file, header=TRUE,sep="\t",row.names="Gene")
colmetadata <- read.table(opt$colmetadata,sep="\t")
rowmetadata <- read.table(opt$rowmetadata,sep="\t")

#metadata_tooltips <- read.table(paste(opt$file,"metadata_tooltips.txt",sep = "."),sep="\t")

#mat <- metadata_tooltips
#mat <- mydata
#mat[] <- paste("\nCountry: ", rowmetadata[, 2], "\nContinent: ", rowmetadata[, 3], "\nCOG: ", colmetadata)

heatmaply(
t(mydata), file = opt$out,
#custom_hovertext = mat,
col_side_colors = colmetadata[, 3],
row_side_colors = rowmetadata[, 2:4],
scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "white" , mid="#E7E4EA", high = "#471777", limits = c(0, 1))
)
