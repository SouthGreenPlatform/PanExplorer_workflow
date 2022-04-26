#!/usr/local/R-4.1.2/bin/R

library(dendextend)
library("optparse")
library(svglite)

#args = commandArgs(trailingOnly=TRUE)

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

svglite(opt$out,width = 31, height = 28)

mydata <- read.table(opt$file, header=TRUE,sep="\t", row.names="Gene")

iris <- mydata

dend_r <- iris %>% dist(method = "man") %>% hclust(method = "ward.D") %>% as.dendrogram %>% ladderize

dend_c <- t(iris) %>% dist(method = "man") %>% hclust(method = "com") %>% as.dendrogram %>% ladderize

mat <- as.matrix(t(iris-1))
out <- gplots::heatmap.2(mat,
          main = "",
          scale="none",
          srtCol=NULL,
          Rowv = dend_c,
          Colv = dend_r,
          key = FALSE,
          margins =c(20,20),
          trace="row", hline = NA, tracecol = NA
         )

write.table(
  data.frame(gene = rownames(mat)[out$rowInd]),
  paste(opt$out, "rows.csv", sep="."),
  row.names = FALSE,
  quote = FALSE,
  sep = ',')

write.table(
  data.frame(gene = colnames(mat)[out$colInd]),
  paste(opt$out, "cols.csv", sep="."),
  row.names = FALSE,
  quote = FALSE,
  sep = ',')

dev.off()

