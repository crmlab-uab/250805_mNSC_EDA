# libraries
library(DESeq2)

## dds
vsd <- vst(dds)
head(colData(vsd))
