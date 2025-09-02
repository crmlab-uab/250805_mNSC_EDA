# libraries
library(DESeq2)
library(BiocParallel)
library(tidyverse)

# Run the main DESeq2 analysis on the protein-coding gene set.
# This single call performs size factor estimation, dispersion estimation, and model fitting.
message("Running DESeq2 analysis on the 'pcg' dataset...")
system.time(des_pcg <- DESeq(dds_list$pcg, parallel = TRUE, BPPARAM = MulticoreParam(nc)))
message("DESeq2 analysis complete.")

# View the coefficient names from the model for constructing contrasts
resultsNames(des_pcg)

# --- Analysis 1: Host effect (NSG vs BL6) within each Driver type ---
message("Extracting results for Host effect within each Driver...")

# Reference level for Driver is "EGFRvIII", so its interaction term is combined with the main effect.
# The main effect "Host_NSG_vs_BL6" is the effect in the reference Driver level.
res_host_by_driver <- list()
drivers <- unique(samples$Driver)

for (drv in drivers) {
  # The main effect is the contrast for the reference level driver
  if (drv == "EGFRvIII") {
    contrast_vec <- c("Host_NSG_vs_BL6")
  } else {
    # For other drivers, we combine the main effect and the interaction term
    contrast_vec <- c("Host_NSG_vs_BL6", paste0("Driver", drv, ".HostNSG"))
  }
  
  res <- results(
    des_pcg,
    contrast = list(contrast_vec),
    alpha = qval,
    lfcThreshold = lfc,
    parallel = TRUE,
    BPPARAM = MulticoreParam(nc)
  )
  
  res_name <- paste0("Host_NSG_vs_BL6_in_", drv)
  res_host_by_driver[[res_name]] <- res
  
  # Summarize and save results
  message(paste("Summary for", res_name))
  summary(res)
  
  # Save RDS object
  saveRDS(res, file = paste0(
    dir_output,
    format(Sys.Date(), "%y%m%d"),
    "_res_",
    res_name,
    ".rds"
  ))
  
  # Save summary CSV
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  write.csv(res_df, file = paste0(
    dir_output,
    format(Sys.Date(), "%y%m%d"),
    "_res_",
    res_name,
    ".csv"
  ))
}

# --- Analysis 2: Driver comparisons within BL6 Host ---
message("Running analysis for Driver comparisons within BL6 host...")

# Subset the dds object for only BL6 samples
dds_bl6 <- dds_list$pcg[, dds_list$pcg$Host == "BL6"]

# Re-level the Driver factor to drop unused levels
dds_bl6$Driver <- droplevels(dds_bl6$Driver)

# Create a new design for this specific comparison and re-run DESeq
design(dds_bl6) <- formula( ~ SeqBatch + Type + Driver)
dds_bl6 <- DESeq(dds_bl6, parallel = TRUE, BPPARAM = MulticoreParam(nc))

res_driver_in_bl6 <- list()
driver_combos <- combn(levels(dds_bl6$Driver), 2) # Get all pairwise combinations

for (i in 1:ncol(driver_combos)) {
  drv1 <- driver_combos[1, i]
  drv2 <- driver_combos[2, i]
  
  res <- results(
    dds_bl6,
    contrast = c("Driver", drv1, drv2),
    alpha = qval,
    lfcThreshold = lfc,
    parallel = TRUE,
    BPPARAM = MulticoreParam(nc)
  )
  
  res_name <- paste0("Driver_", drv1, "_vs_", drv2, "_in_BL6")
  res_driver_in_bl6[[res_name]] <- res
  
  # Summarize and save results
  message(paste("Summary for", res_name))
  summary(res)
  
  # Save RDS object
  saveRDS(res, file = paste0(
    dir_output,
    format(Sys.Date(), "%y%m%d"),
    "_res_",
    res_name,
    ".rds"
  ))
  
  # Save summary CSV
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  write.csv(res_df, file = paste0(
    dir_output,
    format(Sys.Date(), "%y%m%d"),
    "_res_",
    res_name,
    ".csv"
  ))
}

message("All DESeq2 analyses and result extractions are complete.")
