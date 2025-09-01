## libraries
library(tximport)
library(DESeq2)

## count files
txi <- tximport(
  files_sf,
  type = "salmon",
  txIn = TRUE,
  txOut = FALSE,
  tx2gene = tx_gene_symbol
)
names(txi)
head(txi$counts)
all(colnames(txi$counts) == rownames(samples))

## DDS
dds_all <- DESeqDataSetFromTximport(txi, samples, design = des_design)
head(colData(dds_all))

## Estimate size factors and dispersions
dds_all <- estimateSizeFactors(dds_all)
dds_all <- estimateDispersions(dds_all)

## low expression filter
smallest_group_size <- 3
idx <- rowSums(counts(dds_all) >= filt) >= smallest_group_size
dds_filt <- dds_all[idx, ]
nrow(dds_filt)
head(colData(dds_filt))

## Estimate factors
dds_filt <- estimateSizeFactors(dds_filt)

## pcg filter
dds_pcg <- dds_filt[row.names(dds_filt) %in% pcg]
nrow(dds_pcg)
head(colData(dds_pcg))
dds_pcg <- estimateSizeFactors(dds_pcg)

## dds list
dds <- list(dds_all, dds_filt, dds_pcg)
dds.names <- c("all", "filt", "pcg") # nolint: object_name_linter.
names(dds) <- dds.names

for (i in seq_along(dds)) {
  saveRDS(dds[[i]],
          file = paste0(
            dir_output,
            format(Sys.Date(), "%y%m%d"),
            "_",
            "dds",
            "_",
            names(dds[i]),
            ".rds"
          ))
}

## Save raw and norm counts using a loop
dds_types <- paste0("dds_", dds.names)
counts <- list()
cts_names <- c()
for (dds_type in dds_types) {
  for (normalized in c(FALSE, TRUE)) {
    count_type <- ifelse(normalized, "norm", "raw")
    count_data <- data.frame(counts(get(dds_type), normalized = normalized))
    count_name <- paste0("counts_", gsub("dds_", "", dds_type), "_", count_type)
    counts[[count_name]] <- count_data
    cts_names <- c(cts_names, count_name)
  }
}
names(counts) <- cts_names

# Ensure output directory exists
if (!dir.exists(dir_output)) {
  dir.create(dir_output, recursive = TRUE)
}
for (i in names(counts)) {
  write.csv(
    counts[[i]],
    file = paste0(dir_output, format(Sys.Date(), "%y%m%d"), "_", i, ".csv"),
    row.names = TRUE
  )
}

## VST
## transform data using variance stablizing transformation
vsd <- list()
for (i in seq_along(dds)) {
  vsd[[i]] <- vst(dds[[i]],
                  blind = FALSE,
                  nsub = 1000,
                  fitType = "parametric")
}
names(vsd) <- dds.names
for (i in names(vsd)) {
  saveRDS(vsd[[i]], file = paste0(
    dir_output,
    format(Sys.Date(), "%y%m%d"),
    "_vsd",
    "_",
    i,
    ".rds"
  ))
}
