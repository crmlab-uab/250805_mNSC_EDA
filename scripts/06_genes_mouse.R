# libraries
library(biomaRt)

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
