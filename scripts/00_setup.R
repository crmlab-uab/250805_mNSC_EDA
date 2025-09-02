# 00_setup.R
# This script checks for and installs all required R packages.
# Run this once inside your Docker container or R session before knitting the Rmd.

message("Checking for required packages...")

# List of required packages

packages_needed <- c(
  "tidyverse", "rmarkdown", "reactable",   # Core & Reporting
  "BiocManager", "DESeq2", "tximport",    # Bioconductor Core
  "biomaRt", "pheatmap", "RColorBrewer",  # Annotation & Plotting
  "ggplotify"                            # Utility
)

# Function to check and install a package
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (pkg %in% c("DESeq2", "tximport", "biomaRt", "BiocParallel")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  } else {
    message(paste("Package already installed:", pkg))
  }
}

# Install BiocManager first if it's not there
install_if_missing("BiocManager")

# Loop through and install all other packages
for (pkg in packages_needed) {
  install_if_missing(pkg)
}

# A final check for the custom RNAseqQC package
if (!requireNamespace("RNAseqQC", quietly = TRUE)) {
  warning(
    "\nCustom package 'RNAseqQC' not found.\n",
    "The 03_QC.R script will skip the plots that depend on it.\n",
    "To install it, you may need to use remotes::install_github() or remotes::install_local()."
  )
} else {
  message("Custom package 'RNAseqQC' is installed.")
}

message("\nPackage check complete. The environment is ready.")
