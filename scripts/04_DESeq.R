# libraries
library(DESeq2)

# DESeq2
dds_res <- list()
for (i in seq_along(dds_list)) {
  system.time(dds_res[[i]] <- DESeq(dds_list[[i]]))
}
names(dds_res) <- names(dds_list)

# Results
res_list <- list()
for (i in seq_along(dds_res)) {
  system.time(
    res_list[[i]] <-
      results(
        dds_res[[i]],
        contrast = c("condition", "KO", "C"),
        alpha = qval,
        lfcThreshold = lfc,
        parallel = TRUE,
        BPPARAM = MulticoreParam(nc),
        pAdjustMethod = "fdr"))
}