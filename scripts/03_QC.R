# libraries
library(RNAseqQC)
library(pheatmap)

# QC metrics
tc <- list()
lc <- list()
gd <- list()

for (i in seq_along(dds)) {
  tc[[i]] <- plot_total_counts(dds[[i]])
  lc[[i]] <- plot_library_complexity(dds[[i]])
  gd[[i]] <- plot_gene_detection(dds[[i]])
}

names(tc) <- dds.names
names(lc) <- dds.names
names(gd) <- dds.names

ggp <- list(tc, lc, gd)
ggp_names <- c(
  "Total Counts",
  "Library Complexity",
  "Gene Detection"
)
names(ggp) <- ggp_names

# save plots as pdf
for (i in seq_along(ggp)) {
  for (j in seq_along(ggp[[i]])) {
    ggsave(
      filename = paste0(dir_graph, names(ggp[i]), "_", names(ggp[[i]])[j],
                        ".pdf"),
      plot = ggp[[i]][[j]]
    )
  }
}

msd <- list()
for (i in seq_along(vsd)) {
  msd[[i]] <- ggplotify::as.ggplot(mean_sd_plot(vsd[[i]]))
  ggsave(
    filename = paste0(dir_graph, "Mean_SD_", dds.names[i], ".pdf"),
    plot = msd[[i]]
  )
}
names(msd) <- dds.names

### dds QC
histo <- list()
for (i in seq_along(dds)){
  histo[[i]] <- counts(dds[[i]], normalized = TRUE)
}
names(histo) <- dds.names
logcts <- list()
for (i in seq_along(histo)){
  logcts[[i]] <- log(histo[[i]][, 1], 10)
}
d <- list()
for (i in seq_along(logcts)){
  d[[i]] <- density(logcts[[i]])
  plot(d[[i]], main = "",
       xlab = "Normalized Read Counts (log10)", ylab = "Density")
  title(names(histo[i]))
}

# Save and print plotDispEsts results
pdf(file = paste0(dir_graph, "Dispersion_Estimates.pdf"))
plotDispEsts(dds_all)
dev.off()
plotDispEsts(dds_all)
