# @knitr libraries
library(BiocParallel)
library(colorRamps)
library(RColorBrewer)
library(biomaRt)

# @knitr deseq2 variables
pval <- 0.05
qval <- 0.05
lfc <- 1

# @knitr low read counts filter cutoff
filt <- 20

# @knitr date and cores
date <- format(Sys.Date(), format = "%Y%m%d")
nc <- BiocParallel::multicoreWorkers()

# @knitr colors
blp <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ryb <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# @knitr seed
set.seed(888)

# @knitr genome
## /data/project/MillerLab/genomes/GRCm39_vM37

## import tx2gene (contains 3 columns)
## Col 1 - Ensembl txID [transcript_id]
## Col 2 - Ensembl geneID [gene_id]
## Col 3 - gene_name

tx_geneID_genename <- read.delim2(
    file = "/data/input/tx2gene.tsv",
    header = TRUE, sep = "\t")
colnames(tx_geneID_genename)

tx_gene_symbol <- tx_geneID_genename[-2] # remove gene_id column
head(tx_gene_symbol)

tx_gene <- tx_geneID_genename[-3] # remove gene_name column
head(tx_gene)

# ID_to_symbols
genes <- tx_gene_symbol$gene_name
length(genes)
mart <-
  useEnsembl(biomart = "ensembl",
             dataset = "mmusculus_gene_ensembl",
             mirror = "useast")
genes_mouse <-
  biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id_version",
      "mgi_id",
      "mgi_symbol",
      "chromosome_name",
      "strand",
      "start_position",
      "end_position",
      "description",
      "gene_biotype"
    ),
    filters = "external_gene_name",
    values = genes,
    mart = mart
  )
write.csv(genes_mouse,
          file = paste0(
            dir_output,
            format(Sys.Date(), "%y%m%d"),
            "_genes_mouse_vM37.csv"
          ))

# @knitr mouse genes
genes_mouse %>%
  group_by(gene_biotype) %>%
  summarize(n = n())

# Display as an interactive table
reactable(
  genes_mouse,
  searchable = TRUE,
  filterable = TRUE,
  resizable = TRUE
)

# @knitr pcg
pcg <-
  subset(genes_mouse, genes_mouse$gene_biotype == "protein_coding")
dim(pcg)
pcg_df <- pcg

## Remove "RIKEN" genes
pcg_df <- pcg_df[!grepl("RIKEN", pcg_df$description), ]
dim(pcg_df)
## Remove "cDNA sequence" genes
pcg_df <- pcg_df[!grepl("cDNA sequence", pcg_df$description), ]
dim(pcg_df)
## Remove "DNA segment" genes
pcg_df <- pcg_df[!grepl("DNA segment", pcg_df$description), ]
dim(pcg_df)
## Remove "predicted gene" genes
pcg_df <- pcg_df[!grepl("predicted gene", pcg_df$description), ]
dim(pcg_df)
## Remove "MT" chromosomes
pcg_df <- subset(pcg_df, pcg_df$chromosome_name != "MT")
dim(pcg_df)
## Include only somatic chromosomes
pcg_df <- subset(pcg_df, pcg_df$chromosome_name %in% c(1:19, "X", "Y"))

dim(pcg_df)
write.csv(pcg_df, file = paste0(dir_output, "genes_mouse_vM37_pcg.csv"))

pcg <- as.vector(pcg_df$mgi_symbol)

# @knitr genesets
# kinases
file_kinases <- "./genesets/201006_composite_kinases.csv"
if (file.exists(file_kinases)) {
  kinases <- read.csv(file_kinases, header = TRUE, fileEncoding = "UTF-8-BOM")
} else {
  stop(paste("File not found:", file_kinases))
}
str(kinases)
head(kinases)
kinase_genes <- kinases$Mouse_symbol