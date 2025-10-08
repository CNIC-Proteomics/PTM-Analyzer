# Function to install missing packages
install_if_missing <- function(package) {
  if (package == 'limma') {
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    if (!require(package, character.only = TRUE)) {
        BiocManager::install("limma")
    }
    library("limma", character.only = TRUE)
  } else if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "http://cran.us.r-project.org")
    library(package, character.only = TRUE)
  }
}

# List of required packages
packages <- c("limma", "dplyr", "logging", "yaml", "optparse")

# Install missing packages
for (pkg in packages) { 
  install_if_missing(pkg) 
}
