# libraries
library(sftp)

# Set variables
## map local project file to /data
## map local Renviron file to /env
## These mapping assume you have mounted /data and /env as persistent volumes
project_dir <- "/data"
env_dir <- "/env"

# remote_dir assumes you have mounted /data as follows on your local drive
local_dir <- "/Users/rmiller/R/250805_mNSC_EDA"
remote_dir <- "/data/project/MillerLab"

# check RCurl version
RCurl::curlVersion()
options(RCurlOptions = list(verbose = TRUE))

# Set working directory
getwd()
setwd(paste0("/", project_dir))
getwd()

# Load environment
readRenviron(paste0(env_dir, "/.Renviron"))
Sys.getenv("SFTP_USER")
Sys.getenv("SFTP_HOST")
Sys.getenv("SFTP_DIR")

# Check if curl or scp is available for SFTP download
check_tool <- function(tool) {
  suppressWarnings(system(paste("which", tool), intern = TRUE)) != ""
}

download_file <- function(remote, local, user, pass, host, keyfile = NULL) {
  if (check_tool("curl")) {
    message("Using curl for SFTP download...")
    cmd <- sprintf("curl -u %s:%s sftp://%s%s -o %s", user, pass, host, remote, local)
  } else if (check_tool("scp") && !is.null(keyfile)) {
    message("Using scp for SFTP download...")
    cmd <- sprintf("scp -i %s %s@%s:%s %s", keyfile, user, host, remote, local)
  } else {
    stop("Neither curl nor scp is available, or keyfile is missing.")
  }
  system(cmd)
}

# establish sftp connection to Cheaha on /data/project/MillerLab
sftp_con <- sftp_connect(
  server = Sys.getenv("SFTP_HOST"),
  folder = Sys.getenv("SFTP_DIR"),
  username = Sys.getenv("SFTP_USER"),
  password = Sys.getenv("SFTP_PASS"),
  protocol = "sftp://",
  port = Sys.getenv("SFTP_PORT"),
  timeout = 30