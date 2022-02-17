
#!/usr/bin/env Rscript
library("XenofilteR")
args <- commandArgs(trailingOnly = TRUE)
sample.list <- data.frame(human = args[1], mouse=args[2])
bp.param <- SnowParam(workers = 4, type = "SOCK")
XenofilteR(sample.list, destination.folder = args[3], bp.param = bp.param)
