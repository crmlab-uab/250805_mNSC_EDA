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