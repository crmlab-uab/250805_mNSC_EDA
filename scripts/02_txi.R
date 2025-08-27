# libraries
library(tximport)

# make dds object from sf and sample files

## genome - /data/project/MillerLab/genomes/GRCm39_vM37

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

tx2gene <- tx2geneID_genename[-3] # remove gene_name column
head(tx2gene)

# import read count files (salmon output files)
## txi
txi <- tximport(
  files_sf,
  type = "salmon",
  txIn = TRUE,
  txOut = FALSE,
  tx2gene = tx2gene_symbol
)

names(txi)
colnames(txi$counts) <- rownames(samples)
all(colnames(txi$counts) == rownames(samples))