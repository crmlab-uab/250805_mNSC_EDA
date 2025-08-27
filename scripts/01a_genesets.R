# file location
file_kinases <- "./genesets/201006_composite_kinases.csv"

# genesets
kinases <-
  read.csv(
    file_kinases,
    header = TRUE,
    fileEncoding = "UTF-8-BOM"
  )

str(kinases)
head(kinases)

kinase_genes <- kinases$Mouse_symbol