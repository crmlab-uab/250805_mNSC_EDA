# libraries
library(tximport)
library(DESeq2)

# count files
txi <- tximport(
  files_sf,
  type = "salmon",
  tx2gene = tx_gene_symbol
)
names(txi)
head(txi$counts)
all(colnames(txi$counts) == rownames(samples))

# Initialize a master list to hold all dds objects for each model
dds_models <- list()
message("Creating DESeqDataSet objects for each model...")
# Loop through each design formula provided by the parent RMD
for (model_name in names(design_list)) {
  message(paste("... processing model:", model_name))
  design_formula <- design_list[[model_name]]
  
  # Create the initial full DESeqDataSet
  dds_all <- DESeqDataSetFromTximport(txi, samples, design = design_formula)
  
  # Pre-filtering: Keep genes with at least 'filt' counts in the smallest group of samples
  smallest_group_size <- 3
  keep <- rowSums(counts(dds_all) >= filt) >= smallest_group_size
  dds_filt <- dds_all[keep, ]
  
  # Filter for protein-coding genes
  dds_pcg <- dds_filt[rownames(dds_filt) %in% pcg, ]
  
  # Create a named list of the dds objects for the current model
  current_model_dds_list <- list(
    all = dds_all,
    filt = dds_filt,
    pcg = dds_pcg
  )
  
  # Estimate size factors for each dataset within the current model
  current_model_dds_list <- lapply(current_model_dds_list, estimateSizeFactors)
  
  # Add the processed list of dds objects to the master list
  dds_models[[model_name]] <- current_model_dds_list
}

# --- Save Data Objects ---
# Save the master list of dds objects
saveRDS(dds_models, file = paste0(dir_output, date, "_dds_models_list.rds"))

# Perform Variance Stabilizing Transformation (VST) for each model's filtered data
vsd_models <- list()
for (model_name in names(dds_models)) {
  # VST is typically run on the filtered data for QC
  vsd_models[[model_name]] <- vst(dds_models[[model_name]][['pcg']], blind = FALSE)
}
saveRDS(vsd_models, file = paste0(dir_output, date, "_vsd_models_list.rds"))

message("Data preparation complete. DDS and VSD objects saved for all models.")