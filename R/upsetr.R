#!/usr/bin/R

library(UpSetR)
#library(svglite)
args = commandArgs(trailingOnly=TRUE)


movies <- read.csv(file=args[1], header = T, sep = "\t")

nb <- as.integer(args[3])+1

size <- 0.5
if (as.integer(nb)>60){
	size <- 0.2
}

sets <- names(movies[2:nb])
#svglite(args[2], width = 4, height = 4)
pdf(args[2])
#upset(movies,nintersects = 20,number.angles = 0, mb.ratio = c(0.35, 0.65),order.by = "freq", sets = sets, point.size = 0.5, line.size = 0.3, mainbar.y.label = "Intersections", sets.x.label = "Nb genes", text.scale = c(0.8, 0.8, 0.8, 0.8, size, 0.5))
upset(movies,nintersects = 20,number.angles = 0, mb.ratio = c(0.35, 0.65), order.by = "freq",keep.order = T, sets = sets, point.size = 1.5, line.size = 0.3, mainbar.y.label = "Intersections", sets.x.label = "Nb genes", text.scale = c(0.8, 0.8, 0.8, 0.8, size, 0.5))
dev.off()
