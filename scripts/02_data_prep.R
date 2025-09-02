# 02_data_prep.R
message("--- Running 02_data_prep.R: Preparing DESeq2 Objects ---")
library(tximport)
library(DESeq2)

# Import salmon counts
txi <- tximport(
  files_sf,
  type = "salmon",
  tx2gene = tx_gene_symbol,
  ignoreAfterBar = TRUE
)

# Create DESeqDataSet objects for each model
dds_models <- list()
for (model_name in names(design_list)) {
  message(paste("... processing model:", model_name))
  design_formula <- design_list[[model_name]]
  dds_all <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)
  
  keep <- rowSums(counts(dds_all) >= filt) >= 3
  dds_filt <- dds_all[keep, ]
  dds_pcg <- dds_filt[rownames(dds_filt) %in% pcg, ]
  
  current_model_dds_list <- list(all = dds_all, filt = dds_filt, pcg = dds_pcg)
  dds_models[[model_name]] <- lapply(current_model_dds_list, estimateSizeFactors)
}

# Perform Variance Stabilizing Transformation
vsd_models <- list()
for (model_name in names(dds_models)) {
  message(paste("... applying VST to model:", model_name))
  vsd_models[[model_name]] <- vst(dds_models[[model_name]]$filt, blind = TRUE)
}
message("Data preparation complete.")