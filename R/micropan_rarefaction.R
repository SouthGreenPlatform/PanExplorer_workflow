#!/usr/bin/R


library(micropan)
library(dendextend)
library("optparse")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-a", "--alpha"), type="character", default="heaps.tsv",
              help="Heaps output file name [default= %default]", metavar="character"),
  make_option(c("-p", "--pdf"), type="character", default="out.pdf",
              help="PDF output file name [default= %default]", metavar="character")
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
if (is.null(opt$alpha)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (out file).\n", call.=FALSE)
}
if (is.null(opt$pdf)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (out file).\n", call.=FALSE)
}


# Loading a pan-matrix in this package 
#data(xmpl.panmat)

mydata <- read.table(opt$file, header=TRUE,sep="\t", row.names="Gene")
tmydata<-t(mydata)
#xmpl.panmat

rarefaction = rarefaction(tmydata, n.perm = 100)
write.table(
  rarefaction,
  opt$out,
  row.names = FALSE,
  quote = FALSE,
  sep = '\t')


# Estimating population openness
h.est <- heaps(tmydata, n.perm = 100)
write.table(
  h.est,
  opt$alpha,
  row.names = TRUE,
  quote = FALSE,
  sep = '\t')

# If alpha < 1 it indicates an open pan-genome

pdf(opt$pdf)

library(ggplot2)
library(tidyr)
rarefaction %>%
gather(key = "Permutation", value = "Clusters", -Genome) %>%
ggplot(aes(x = Genome, y = Clusters, group = Permutation)) +
geom_point() +
stat_summary(aes(y = Clusters,group=1), fun=mean, colour="red", geom="line",group=1)
