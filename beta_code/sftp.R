# libraries
library(sftp)

# Set variables
## map local project file to /data
## map local Renviron file to /env
## These mapping assume you have mounted /data and /env as persistent volumes
project_dir <- "/data" # mapped to local R folder

# remote_dir assumes you have mounted /data as follows on your local drive
local_dir <- paste0(project_dir, "/250805_mNSC_EDA")

# check RCurl version
RCurl::curlVersion()
options(RCurlOptions = list(verbose = TRUE))

# confirm wd is /home/rstudio
getwd()

# Load environment
## This code assumes you have placed your .Renviron file in the root directory 
## that is mapped to /data
readRenviron(paste0(getwd(), "/.Renviron"))
Sys.getenv("SFTP_USER")
Sys.getenv("SFTP_HOST")
Sys.getenv("SFTP_DIR")

# Set working directory
setwd(project_dir)
getwd()

# establish sftp connection to Cheaha on /data/project/MillerLab
sftp_con <- sftp_connect(
  server = Sys.getenv("SFTP_HOST"),
  folder = Sys.getenv("SFTP_DIR"),
  username = Sys.getenv("SFTP_USER"),
  protocol = "sftp://",
  port = Sys.getenv("SFTP_PORT"),
  timeout = 60
)

# Run download
RCurl::listCurlOptions()

yo <- sftp_list(
  sftp_connection = sftp_con,
  verbose = TRUE,
  curlPerformVerbose = TRUE,
  recurse = TRUE,
  curl_options = list(
    ftp.ssl = TRUE,
    ssh.private.keyfile = Sys.getenv("SFTP_PRIVATE_KEY"),
    ssh.public.keyfile = Sys.getenv("SFTP_PUBLIC_KEY")
  ))

file_kinases <- "./genesets/201006_composite_kinases.csv"

## genesets
kinases <-
  read.csv(
    file_kinases,
    header = TRUE,
    fileEncoding = "UTF-8-BOM"
  )

str(kinases)
head(kinases)

kinase_genes <- kinases$HGNC