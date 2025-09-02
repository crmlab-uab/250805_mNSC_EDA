# 03_QC.R
message("--- Running 03_QC.R: Quality Control ---")
library(pheatmap)
library(ggplot2)

# --- RNAseqQC plots (will run only if the package is installed) ---
if (requireNamespace("RNAseqQC", quietly = TRUE)) {
  message("RNAseqQC package found. Generating standard QC plots...")
  qc_plots <- list()
  # Corrected to loop through dds_models
  for (name in names(dds_models)) {
    qc_plots[[name]][['total_counts']] <- RNAseqQC::plot_total_counts(dds_models[[name]]$filt)
    qc_plots[[name]][['library_complexity']] <- RNAseqQC::plot_library_complexity(dds_models[[name]]$filt)
    qc_plots[[name]][['gene_detection']] <- RNAseqQC::plot_gene_detection(dds_models[[name]]$filt)
  }
  # Save plots to PDF
  for (name in names(qc_plots)) {
    for (plot_type in names(qc_plots[[name]])) {
      ggsave(
        filename = paste0(dir_graph, plot_type, "_", name, ".pdf"),
        plot = qc_plots[[name]][[plot_type]]
      )
    }
  }
} else {
  message("RNAseqQC package not found. Skipping plots that depend on it.")
}

# --- Mean-SD plots ---
for (name in names(vsd_models)) {
  pdf(paste0(
    dir_graph,
    format(Sys.Date(), "%y%m%d"),
    "_Mean_SD_",
    name,
    ".pdf"
  ))
  DESeq2::meanSdPlot(assay(vsd_models[[name]]))
  dev.off()
}

# Density plots of normalized counts
for (name in names(dds_models)) {
  norm_counts <- counts(dds_models[[name]], normalized = TRUE)
  # Adding +1 to avoid log(0).
  log_norm_counts <- log10(norm_counts + 1)
  # Convert to a long format for ggplot
  df_long <- as.data.frame(log_norm_counts) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "sample",
                        values_to = "log10_counts")
  p <- ggplot(df_long, aes(x = log10_counts, color = sample)) +
    geom_density() +
    labs(title = paste("Count Density:", name),
         x = "Normalized Read Counts (log10)",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "none") # Remove legend if too crowded
  ggsave(filename = paste0(dir_graph, "DensityPlot_", name, ".pdf"),
         plot = p)
}

# Dispersion Estimates Plot (run only on the object used for DEA)
pdf(file = paste0(dir_graph, "Dispersion_Estimates.pdf"))
plotDispEsts(dds_models$pcg)
dev.off()

# --- Boxplots of Count Distributions ---
# This QC step is crucial for visualizing library size differences and the effect of normalization.
message("Generating QC Boxplots of counts...")
for (model_name in names(dds_models)) {
  dds <- dds_models[[model_name]][['pcg']] # Use the protein-coding filtered data
  
  # Prepare data for plotting
  raw_counts <- counts(dds, normalized = FALSE)
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Combine into a single data frame for ggplot
  df_raw <- as.data.frame(log2(raw_counts + 1)) %>%
    tidyr::pivot_longer(everything(), names_to = "sample", values_to = "log2_count") %>%
    mutate(type = "Raw")
  
  df_norm <- as.data.frame(log2(norm_counts + 1)) %>%
    tidyr::pivot_longer(everything(), names_to = "sample", values_to = "log2_count") %>%
    mutate(type = "Normalized")
  
  df_counts <- rbind(df_raw, df_norm)
  df_counts$type <- factor(df_counts$type, levels = c("Raw", "Normalized"))
  
  p <- ggplot(df_counts, aes(x = sample, y = log2_count, fill = type)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap( ~ type, scales = "free_y", ncol = 1) +
    labs(
      title = paste("Count Distributions (Model:", model_name, ")"),
      x = "Sample",
      y = "Log2(count + 1)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))
  
  ggsave(
    filename = paste0(dir_graph, date, "_QC_Boxplot_Counts_", model_name, ".pdf"),
    plot = p,
    width = 12,
    height = 8
  )
}

# --- PCA Plots ---
message("Generating QC PCA plots...")
for (model_name in names(vsd_models)) {
  pca_plot <- plotPCA(vsd_models[[model_name]], intgroup = c("Driver", "Host", "Type")) +
    labs(title = paste("PCA on VST (Model:", model_name, ")")) + theme_bw()
  print(pca_plot)
}

# --- Sample-to-Sample Distance Heatmaps ---
message("Generating QC Sample Distance Heatmaps...")
for (model_name in names(vsd_models)) {
  vsd <- vsd_models[[model_name]]
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colors,
    annotation_col = anno,
    annotation_colors = anno_colors,
    main = paste("Sample-to-Sample Distance (Model:", model_name, ")"),
    filename = paste0(dir_graph, date, "_QC_SampleDistanceHeatmap_", model_name, ".pdf"),
    width = 10, height = 8
  )
}

message("QC script finished.")