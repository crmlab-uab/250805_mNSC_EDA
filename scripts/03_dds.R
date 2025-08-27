# libraries
library(DESeq2)

## dds
dds <- DESeqDataSetFromTximport(
  txi,
  samples,
  design = ~ Type + Driver + Host)

head(colData(dds))

# Estimate size factors
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
