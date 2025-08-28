# libraries
library(tximport)
library(DESeq2)

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
head(txi$counts)
all(colnames(txi$counts) == rownames(samples))

## dds
dds_all <- DESeqDataSetFromTximport(txi, samples, design = ~ Type + Driver + Host)
head(colData(dds_all))

## Estimate size factors and dispersions
dds_all <- estimateSizeFactors(dds_all)
dds_all <- estimateDispersions(dds_all)

# low expression filter
keep <- rowSums(counts(dds_all)) > 5
dds_filt <- dds_all[keep, ]
nrow(dds_filt)
head(colData(dds_filt))

## Estimate factors
dds_filt <- estimateSizeFactors(dds_filt)
normalizationFactors(dds_filt)

# pcg filt
dds_pcg <- dds_filt[row.names(dds_filt) %in% pcg]
nrow(dds_pcg)
head(colData(dds_pcg))
dds_pcg <- estimateSizeFactors(dds_pcg)
normalizationFactors(dds_pcg)
saveRDS(dds_pcg, file = paste0(dir_output, format(Sys.Date(), "%y%m%d"), "_", "dds_pcg", ".rds"))

# dds list
dds <- list(dds_all, dds_filt, dds_pcg)
dds_names <- c("dds_all", "dds_filt", "dds_pcg")
names(dds) <- dds_names

#Save raw and norm counts
counts_all_norm <- data.frame(counts(dds_all, normalized = TRUE))
counts_all_raw <- data.frame(counts(dds_all, normalized = FALSE))
counts_filt_norm <- data.frame(counts(dds_filt, normalized = TRUE))
counts_filt_raw <- data.frame(counts(dds_filt, normalized = FALSE))
counts_pcg_norm <- data.frame(counts(dds_pcg, normalized = TRUE))
counts_pcg_raw <- data.frame(counts(dds_pcg, normalized = FALSE))

counts <- list(
  counts_all_raw,
  counts_all_norm,
  counts_filt_raw,
  counts_filt_norm,
  counts_pcg_raw,
  counts_pcg_norm
)
cts_names <- c(
  "counts_all_raw",
  "counts_all_norm",
  "counts_filt_raw",
  "counts_filt_norm",
  "counts_pcg_raw",
  "counts_pcg_norm"
)
names(counts) <- cts_names

for (i in seq_along(counts)) {
  write.csv(counts[[i]],
            file = paste0(dir_output, names(counts[i]), ".csv"),
            row.names = TRUE)
}

# VST
## transform data using vst
vsd <- list()
for (i in seq_along(dds)){
  vsd[[i]] <- vst(dds[[i]],
                  blind = FALSE,
                  nsub = 1000,
                  fitType = "parametric")
}
vsd_names <- c("vsd_all", "vsd_filt", "vsd_pcg")
names(vsd) <- vsd_names

# remove counts, dds, and vsd objects
objects_to_remove <- ls(pattern = "counts_")
rm(list = objects_to_remove)

objects_to_remove <- ls(pattern = "dds_")
rm(list = objects_to_remove)

objects_to_remove <- ls(pattern = "vsd_")
rm(list = objects_to_remove)