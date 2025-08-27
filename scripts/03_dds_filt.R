# libraries
library(DESeq2)

# low expression filter
keep <- rowSums(counts(dds)) > 5
dds_filt <- dds[keep,]
nrow(dds_filt)
head(colData(dds_filt))

# Estimate size factors
dds_filt <- estimateSizeFactors(dds_filt)

normalizationFactors(dds_filt)

# pcg filt
## filter for PCG only
dds_pcg <- dds_filt[row.names(dds_filt) %in% pcg]
nrow(dds_pcg)
head(colData(dds_pcg))

# Estimate size factors
dds_pcg <- estimateSizeFactors(dds_pcg)

normalizationFactors(dds_pcg)

# Save dds_ as RDS
saveRDS(dds_pcg, file = paste0(dir_output,format(Sys.Date(),"%y%m%d"),"_","dds_pcg",".rds"))

#Save raw and norm counts after filtering
counts_filt_norm <- as.data.frame(counts(dds_filt, normalized=TRUE))
counts_filt_raw <- as.data.frame(counts(dds_filt, normalized=FALSE))

counts_pcg_norm <- as.data.frame(counts(dds_pcg, normalized=TRUE))
counts_pcg_raw <- as.data.frame(counts(dds_pcg, normalized=FALSE))

write.csv(
  counts_filt_norm,
  file=paste0(dir_output,"counts_filt_norm",".csv"),
  row.names = TRUE
)

write.csv(
  counts_filt_raw,
  file=paste0(dir_output,"counts_filt_raw",".csv"),
  row.names = TRUE
)

write.csv(
  counts_pcg_norm,
  file=paste0(dir_output,"counts_pcg_norm",".csv"),
  row.names = TRUE
)

write.csv(
  counts_pcg_raw,
  file=paste0(dir_output,"counts_pcg_raw",".csv"),
  row.names = TRUE
)

# VST
# transform filtered data using vst
vsd_filt<-vst(dds_filt, blind=FALSE, nsub=1000, fitType = "parametric")
