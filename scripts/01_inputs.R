# 01_inputs.R
message("--- Running 01_inputs.R: Loading Data ---")
library(biomaRt)
library(reactable)

# Load sample metadata from the project root directory
samples <- read.csv("./input/samples.csv", row.names = 1, header = TRUE)
message("Displaying sample metadata:")
print(reactable(samples))

# Load tx2gene mapping file
tx2gene_file <- "./input/tx2gene.tsv"
if (!file.exists(tx2gene_file)) {
  stop("tx2gene mapping file not found at: ", tx2gene_file)
}
tx_geneID_genename <- read.delim2(file = tx2gene_file, header = TRUE, sep = "\t")
tx_gene_symbol <- tx_geneID_genename[, c(1, 3)]
colnames(tx_gene_symbol) <- c("TXNAME", "GENEID")

# Get paths to salmon quantification files, assuming they are in sub-folders
files_sf <- paste0("./input/", rownames(samples), ".sf")
names(files_sf) <- rownames(samples)
if (!all(file.exists(files_sf))) {
  stop("One or more Salmon quantification files are missing.")
}
message("Salmon quantification file paths generated.")

# Gene Annotation using biomaRt
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes <-
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

write.csv(genes,
          file = paste0(
            dir_output,
            format(Sys.Date(), "%y%m%d"),
            "_genes_mouse_vM37.csv"
          ))

# mouse genes
genes %>%
  group_by(gene_biotype) %>%
  summarize(n = n())

# pcg
pcg_df <-
  subset(genes_mouse, genes_mouse$gene_biotype == "protein_coding")
dim(pcg_df)

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
message("Gene annotations fetched.")