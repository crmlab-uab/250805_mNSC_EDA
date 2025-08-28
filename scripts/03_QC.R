# libraries
library(RNAseqQC)
library(DESeq2)
library(pheatmap)

# QC metrics
tc <- list()
lc <- list()
gd <- list()
disp <- list()

for (i in seq_along(dds)) {
  tc[[i]] <- plot_total_counts(dds[[i]])
  lc[[i]] <- plot_library_complexity(dds[[i]])
  gd[[i]] <- plot_gene_detection(dds[[i]])
  disp[[i]] <- plotDispEsts(dds[[i]])
}

ggp <- list(tc, lc, gd, disp)
ggp_names <- c("Total Counts", "Library Complexity", "Gene Detection", "Dispersion Estimates")
names(ggp) <- ggp_names

for (i in 1:length(dds)) {
  for (j in 1:length(ggp)) {
    names(ggp[i][j]) <- names(dds[i])
  }
}

msd <- list()
for (i in seq_along(vsd)){
  msd[[i]] <- mean_sd_plot(vsd[[i]])
  walk(msd, print)
}
msd_names <- c("all", "filtered", "pcg")
names(msd) <- msd_names

######### START HERE #######


### dds_filter QC
histo <- counts(dds_filt, normalized=TRUE)
logcounts <- log(histo[,1],10) 
d <- density(logcounts)
plot(d,main="",xlab="Normalized Read Counts (log10)", ylab="Density")

for (s in 1:nrow(samples)){
  logcounts <- log(histo[,s],10) 
  d <- density(logcounts)
  lines(d)
}

logcounts <- log(histo[,1],10) 
d <- density(logcounts)
plot(d,main="",xlab="Normalized Read Counts (log10)", ylab="Density")

for (s in 1:nrow(samples)){
  logcounts <- log(histo[,s],10) 
  d <- density(logcounts)
  lines(d)
}

### dds_pcg QC
histo <- counts(dds_pcg, normalized=TRUE)
logcounts <- log(histo[,1],10) 
d <- density(logcounts)
plot(d,main="",xlab="Normalized Read Counts (log10)", ylab="Density")

for (s in 1:nrow(samples)){
  logcounts <- log(histo[,s],10) 
  d <- density(logcounts)
  lines(d)
}

logcounts <- log(histo[,1],10) 
d <- density(logcounts)
plot(d,main="",xlab="Normalized Read Counts (log10)", ylab="Density")

for (s in 1:nrow(samples)){
  logcounts <- log(histo[,s],10) 
  d <- density(logcounts)
  lines(d)
}

# sample euclidian distances filtered data
sampleDists <- dist(t(assay(vsd_filt)))
sampleDistMatrix <- as.matrix(sampleDists)

# plot sample euclidian distances via a heatmap
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = anno,
  fontsize = 8,
  show_colnames = FALSE,
  scale = "row",
  show_rownames = FALSE,
  color = blp
)

# save to PDF
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
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
pcaData <- plotPCA(vsd, intgroup = c(var1, var2, var3), returnData = TRUE)
write.csv(pcaData,
          paste0(dir_output,
                 "pcaData.csv"
          ))
percentVar <- round(100 * attr(pcaData, "percentVar"))

for (i in c(var1, var2, var3)){
p[[i]] <- ggplot(pcaData, aes(x = PC1, y = PC2, color = i)) +
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