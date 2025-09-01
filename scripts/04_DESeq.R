## libraries
library(DESeq2)
library(BiocParallel)

## DESeq2
system.time(des_pcg <- DESeq(dds$pcg, parallel = TRUE, BPPARAM = MulticoreParam(nc)))

des_pcg_bl6 <- des_pcg

## results
system.time(
  res <-
    results(
      des_pcg,
      contrast = c("Host", "BL6", "NSG"),
      alpha = qval,
      lfcThreshold = lfc,
      parallel = TRUE,
      BPPARAM = MulticoreParam(nc),
      pAdjustMethod = "fdr"
    )
)
summary(res)