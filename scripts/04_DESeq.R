# libraries
library(DESeq2)
library(BiocParallel)

# DESeq2
system.time(des_pcg <- DESeq(dds$pcg, parallel = TRUE, BPPARAM = MulticoreParam(nc)))

# DES subsets
des_pgc_subset_1 <- subset(des_pcg, rownames(des_pcg) %in% rownames(subset_1))
des_pgc_subset_2 <- subset(des_pcg, rownames(des_pcg) %in% rownames(subset_2))

# results
res_subset_1 <- list()
system.time(for (i in unique(subset_1$Driver)) {
  res_subset_1[[i]] <-
    results(
      des_pgc_subset_1,
      contrast = c("Host", "BL6", "NSG"),
      alpha = qval,
      lfcThreshold = lfc,
      parallel = TRUE,
      BPPARAM = MulticoreParam(nc),
      pAdjustMethod = "fdr"
    )
})
summary(res)