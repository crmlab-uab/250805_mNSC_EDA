# libraries
library(sftp)
library(RCurl)
library(tximport)
library(DESeq2)
library(RNAseqQC)

# make dds object from sf and sample files

# import tx2gene (contains 3 columns)
  # Ensembl txID [transcript_id]
  # Ensembl geneID [gene_id]
  # gene_name

tx2geneID_genename <- read.delim2(
  file = "/data/input/tx2gene.tsv", 
  header = TRUE, sep = '\t')
colnames(tx2geneID_genename)

tx2gene_symbol <- tx2geneID_genename[-2] # remove gene_id column
head(tx2gene_symbol)

# import read count files (salmon output files)

sftp_con <- sftp_connect(
  server = "cheaha.rc.uab.edu",
  folder = "/data/project/MillerLab/projects/mNSC/250730_mNSC_nf-core/",
  username = "rmiller",
  password = "T3nn3ss33V0l$1998#1",
  protocol = "sftp://",
  port = 22,
  timeout = 300
)


## txi

txi <- tximport(files_sf, type = "salmon", tx2gene = tx2gene_symbol)

names(txi)
colnames(txi$counts) <- rownames(samples)
all(colnames(txi$counts) == rownames(samples))

counts_raw <- as.data.frame(txi$counts)
look_for <- c("Egfr", "Pdgfra", "Pten", "Cdkn2a")
counts_raw[rownames(counts_raw) %in% look_for, ]

## dds
dds <- DESeqDataSetFromTximport(
  txi,
  samples,
  design = ~ var7, var6, var5, var4, var3, var2, var1)

# Estimate size factors
dds <- estimateSizeFactors(dds)

# RNAseqQC metrics
plot_total_counts(dds)
plot_library_complexity(dds)
plot_gene_detection(dds)
