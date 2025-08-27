# libraries
library(BiocParallel)
library(colorRamps)
library(RColorBrewer)

## variables
pval <- 0.05
qval <- 0.05
lfc <- 1

date <- format(Sys.Date(), format = "%Y%m%d")
nc <- BiocParallel::multicoreWorkers()

blp <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ryb <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

## folders

# output
dir_output <- "./output/"
dir.create(paste0(dir_output))

# graphs
dir_graph <- "./graphs/"
dir.create(paste0(dir_graph))

## seed
set.seed(888)