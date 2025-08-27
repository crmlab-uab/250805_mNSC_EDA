# libraries
library(RNAseqQC)
library(DESeq2)

# RNAseqQC metrics
p <- plot_total_counts(dds)
print(p)
p <- plot_library_complexity(dds)
print(p)
p <- plot_gene_detection(dds)
print(p)
p <- plotDispEsts(dds)
print(p)
p <- mean_sd_plot(vsd)
print(p)

# RNAseqQC metrics - filt
p <- plot_total_counts(dds_filt)
print(p)
p <- plot_library_complexity(dds_filt)
print(p)
p <- plot_gene_detection(dds_filt)
print(p)
p <- plotDispEsts(dds_filt)
print(p)
p <- mean_sd_plot(vsd_filt)
print(p)

# RNAseqQC metrics - pcg
p <- plot_total_counts(dds_pcg)
print(p)
p <- plot_library_complexity(dds_pcg)
print(p)
p <- plot_gene_detection(dds_pcg)
print(p)
p <- plotDispEsts(dds_pcg)
print(p)