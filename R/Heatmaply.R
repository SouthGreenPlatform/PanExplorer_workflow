#!/usr/bin/R

library("optparse")
library(heatmaply)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
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

mydata <- read.table(opt$file,sep="\t",fill=TRUE,header=TRUE, row.names = 1)
heatmaply(mydata,file = "heatmaply.html",plot_method="plotly",scale_fill_gradient_fun = ggplot2::scale_fill_gradient2( low = "white" , high = "blue", limits = c(0, 100)))
