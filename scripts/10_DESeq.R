# libraries
library(DESeq2)

des <- DESeq(dds_filt)

# QC
plotPCA(vsd_filt,  intgroup = c("Host", "Driver", "Type"))