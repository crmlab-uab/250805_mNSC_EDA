# libraries
library(RNAseqQC)
library(DESeq2)

# RNAseqQC metrics
p <- plot_total_counts(dds_pcg)
print(p)
p <- plot_library_complexity(dds_pcg)
print(p)
p <- plot_gene_detection(dds_pcg)
print(p)
p <- plotDispEsts(dds_pcg)
print(p)