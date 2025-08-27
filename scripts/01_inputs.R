# libraries
library(BiocParallel)
library(colorRamps)
library(RColorBrewer)
library(biomaRt)

# variables
pval <- 0.05
qval <- 0.05
lfc <- 1

date <- format(Sys.Date(), format = "%Y%m%d")
nc <- BiocParallel::multicoreWorkers()

blp <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ryb <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# folders
dir_output <- "./output/"
dir.create(paste0(dir_output))

dir_graph <- "./graphs/"
dir.create(paste0(dir_graph))

# seed
set.seed(888)

## genome - /data/project/MillerLab/genomes/GRCm39_vM37

# import tx2gene (contains 3 columns)
# Col 1 - Ensembl txID [transcript_id]
# Col 2 - Ensembl geneID [gene_id]
# Col 3 - gene_name

tx2geneID_genename <- read.delim2(
  file = "/data/input/tx2gene.tsv", 
  header = TRUE, sep = '\t')
colnames(tx2geneID_genename)

tx2gene_symbol <- tx2geneID_genename[-2] # remove gene_id column
head(tx2gene_symbol)

tx2gene <- tx2geneID_genename[-3] # remove gene_name column
head(tx2gene)

## ID2Symbols
genes <- tx2gene_symbol$gene_name
length(genes)

mart <-
  useEnsembl(biomart = "ensembl",
             dataset = "mmusculus_gene_ensembl",
             mirror = "useast")

genes_mouse <-
  biomaRt::getBM(
    attributes = c("ensembl_gene_id_version",
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
    mart=mart
  )
write.csv(genes_mouse, file = paste0(dir_output, "genes_mouse_vM37.csv"))

genes_mouse %>%
  group_by(gene_biotype) %>%
  summarize(n=n())

## pcg
pcg <-
  subset(genes_mouse,
         genes_mouse$gene_biotype == "protein_coding")
dim(pcg)
pcg_df <- pcg

# Remove "RIKEN" genes by partial string match
pcg_df <- pcg_df[!grepl("RIKEN",pcg_df$description),]
dim(pcg_df)
# Remove "cDNA sequence" genes by partial string match
pcg_df <- pcg_df[!grepl("cDNA sequence",pcg_df$description),]
dim(pcg_df)
# Remove "DNA segment" genes by partial string match
pcg_df <- pcg_df[!grepl("DNA segment",pcg_df$description),]
dim(pcg_df)
# Remove "predicted gene" genes by partial string match
pcg_df <- pcg_df[!grepl("predicted gene",pcg_df$description),]
dim(pcg_df)
# Remove "MT" chromosomes
pcg_df <- subset(pcg_df, pcg_df$chromosome_name!="MT")
dim(pcg_df)
# Include only chromosomes
pcg_df <- subset(pcg_df, pcg_df$chromosome_name %in% c(1:19, "X", "Y"))
dim(pcg_df)

pcg <- as.vector(pcg_df$mgi_symbol)

write.csv(pcg_df, file = paste0(dir_output, "genes_mouse_vM37_pcg.csv"))

# import genesets
## kinases
file_kinases <- "./genesets/201006_composite_kinases.csv"
kinases <-
  read.csv(file_kinases, header = TRUE, fileEncoding = "UTF-8-BOM")
str(kinases)
head(kinases)
kinase_genes <- kinases$Mouse_symbol