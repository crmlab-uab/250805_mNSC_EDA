# libraries
library(tximport)
library(DESeq2)

# make dds object from sf and sample files

# import read count files (salmon output files)
## txi
txi <- tximport(
  files_sf,
  type = "salmon",
  txIn = TRUE,
  txOut = FALSE,
  tx2gene = tx2gene_symbol
)
names(txi)
head(txi$counts)
all(colnames(txi$counts) == rownames(samples))

## dds
dds <- DESeqDataSetFromTximport(
  txi,
  samples,
  design = ~ Type + Driver + Host)
head(colData(dds))

# Estimate size factors and dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

## vst
vsd <- vst(dds)
head(colData(vsd))