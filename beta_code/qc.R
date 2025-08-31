
# sample euclidian distances filtered data
sample_dist <- dist(t(assay(vsd_filt)))
sample_dist_matrix <- as.matrix(sample_dist)

# plot sample euclidian distances via a heatmap
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  annotation_col = anno,
  fontsize = 8,
  show_colnames = FALSE,
  scale = "row",
  show_rownames = FALSE,
  color = blp
)

# save to PDF
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  annotation_col = anno,
  annotation_colors = anno_colors,
  fontsize = 8,
  show_colnames = FALSE,
  scale = "row",
  show_rownames = FALSE,
  color = blp,
  filename = paste(dir_graph, "sampleDistEucl",
                   "pdf",
                   sep = ".")
)

# PCA plots (vsd_filt)
pca_data <- plotPCA(vsd, intgroup = c(var1, var2, var3), returnData = TRUE)
write.csv(pca_data,
          paste0(dir_output,
                 "pca_data.csv" # nolint: indentation_linter.
          ))
percentVar <- round(100 * attr(pca_data, "percentVar"))

for (i in c(var1, var2, var3)){
  p[[i]] <- ggplot(pca_data, aes(x = PC1, y = PC2, color = i)) +
    geom_point(size =3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed()
  
  p[[i]]
  # save graph as PDF
  ggsave(
    filename = paste0(dir_graph, "pca_vsd_filt", i, ".pdf"),
    plot =
      p[[i]]
  )
}

p <- PCAtools::pca(vsd_filt, metadata = samples, removeVar = 0.1)
screeplot(p)

head(p$loadings)
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 2.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3',
             caption = 'Top 1% variables',
             shape = 24,
             colMidpoint = 0,
             col = c("limegreen", "black", "red3"),
             drawConnectors = TRUE,
             returnPlot = TRUE)
ggsave(filename = paste0(dir_graph, "loadings_plot_vsd_filt.pdf"), plot = last_plot())

# RNAseqQC metrics
pca_res <- plot_pca(vsd_filt, show_plot = FALSE)
plot_loadings(pca_res, PC = 1, annotate_top_n = 6)
plot_loadings(pca_res, PC = 1, highlight_genes = project_genes)

# loadings
rv = rowVars(assay(vsd_filt)) 
select = order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc = prcomp(t(assay(vsd_filt)[select,]))
loadings = as.data.frame(pc$rotation)
aload = abs(loadings)
sweep(aload, 2, colSums(aload), "/")
View(aload)