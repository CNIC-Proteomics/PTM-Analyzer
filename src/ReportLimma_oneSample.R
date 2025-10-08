# more information
'
https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
16.1 Swirl Zebrafish: A Single-Group Experiment
pg. 76
'

# ignore warnings
suppressWarnings({
  library(yaml)
})

# import modules
library("limma")
library("dplyr")
library("logging")
library("yaml")
library("optparse")

# define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path to the isanxot report file", metavar = "character"),
  make_option(c("-s", "--samples"), type = "character", default = NULL, help = "Path to the sample table file", metavar = "character"),
  make_option(c("-c", "--config"), type = "character", default = NULL, help = "Path to the YAML configuration file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Path to the output file", metavar = "character")
)

# parse options
opt_parser <- OptionParser(
  option_list = option_list,
  description = paste(
    "ReportLimma_oneSample: Performs the analysis of differential expression using Limma for a one-sample test.",
    sep = "\n"
  )
)
opt <- parse_args(opt_parser)

# ensure required arguments are provided
if (is.null(opt$config) || is.null(opt$input) || is.null(opt$samples) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments!", call. = FALSE)
}
# opt <- list()
# opt$input    <- "S:/U_Proteomica/LABS/LAB_DF/PTM_Analysis_Comet/FA_Biopsias/2_ReportLimma_wo_GUI_oneSample_DEL/NM_Tabla_final.tsv"
# opt$samples  <- "S:/U_Proteomica/LABS/LAB_DF/PTM_Analysis_Comet/FA_Biopsias/2_ReportLimma_wo_GUI_oneSample_DEL/limma_comparisons.tsv"
# opt$config   <- "S:/U_Proteomica/LABS/LAB_DF/PTM_Analysis_Comet/FA_Biopsias/2_ReportLimma_wo_GUI_oneSample_DEL/params.yml"
# opt$output   <- "S:/U_Proteomica/LABS/LAB_DF/PTM_Analysis_Comet/FA_Biopsias/2_ReportLimma_wo_GUI_oneSample_DEL/LIMMA_NM_Tabla_final.tsv"




calculatePvalues <- function(headers, data, config, sampleGroups) {
  
  # init output
  # list of tuple [(low.LEVEL, df_data, df_headers),()...]
  output <- list()
  
  index_output <- 1
  index_intg <- 1
  for (integration in config$integrations) {

    # filter data based on the current integration
    integration_idx <- which(headers[1, ] %in% integration) # extract indices of columns in which row 1 matches any target
    reportData_intg <- data[,integration_idx]
    # get the headers for the integration
    headers_intg <- headers[,integration_idx, drop = FALSE]
    colnames(reportData_intg) <- headers_intg[2,]
    
    # filter data based on lowLevel
    # get vector with low level
    lowLevel1 <- config[["lowLevel1"]][index_intg]
    lowLevel2 <- config[["lowLevel2"]][index_intg]
    lowLevel_idx <- which(headers[1, ] %in% lowLevel1 & headers[2, ] %in% lowLevel2)
    reportData_low <- data[, lowLevel_idx, drop = FALSE]
    # get the headers for the lowLevel
    headers_low <- headers[, lowLevel_idx, drop = FALSE]
    colnames(reportData_low) <- headers_low[1,]

    # combine data with the lowLevel and the data with the integration
    reportData2 <- cbind(reportData_low, reportData_intg)
    headers2 <- cbind(headers_low, headers_intg)

    # go through the groups
    for (group in names(sampleGroups)) {
      # get the samples of current group
      samples <- sampleGroups[[group]]
      loginfo(paste0("Prepare dataset for '",integration,"' > '", group, "': [", paste(samples, collapse = ","), "]"), logger="ReportLimma")

      # create design for the group
      design <- matrix(rep(1, length(samples)))
      rownames(design) <- samples
      colnames(design) <- group
      
      # get the data for the current LowLevel and the samples
      df <- reportData2[, c(lowLevel1, samples)]
      df <- df[!duplicated(df), ] # remove duplicated rows
      rownames(df) <- df[,lowLevel1] # indexed the lowLevel values
      df[,lowLevel1] <- NULL
      
      # calculate Limma (one sample)
      loginfo(paste0("Calculating limma for '", lowLevel1, "' > '", group, "'..."), logger="ReportLimma")
      fit <- lmFit(df,design)
      fit <- eBayes(fit)
      pvalue <- as.data.frame(fit$p.value)
      meanRow <- as.data.frame(fit$Amean)
      limma_output <- merge(meanRow, pvalue, 'row.names', all=T)
      colnames(limma_output) <- c(lowLevel1, 'dX', 'pvalue')
      
      # calculate LPS
      limma_output$LPS <- -log10(limma_output$pvalue) * sign(limma_output$dX)
      
      # rename with integration prefix (except lowLevel)
      limma_header_row1 <- c(
                    paste0(integration, "_dX_", group), 
                    paste0(integration, "_limma_", group),
                    paste0(integration, "_logLimma_", group))
      limma_header_row2 <- c("dX", "pvalue", "LPS")
      limma_headers <- rbind(limma_header_row1, limma_header_row2)
      
      # append limma outputs into global output
      output[[index_output]] <- list(paste0(lowLevel1,".",lowLevel2), limma_output, limma_headers)
      
      # increase index of ouptut (intg/sample) loop
      index_output <- index_output +1
    }
    
    # increase index of integration loop
    index_intg <- index_intg +1
  }
  
  return (output)
  
}

mergeOutputList <- function (outputList, headers, data) {
  
  # get output and output header
  output_headers <- headers
  output <- data
  # create unique header for output
  colnames(output) <- paste0(headers[1, ], ".", headers[2, ])
  
  for (i in 1:length(outputList)) {
    # get limma outputs
    limma_lowLevel <- outputList[[i]][[1]]
    limma_output   <- outputList[[i]][[2]]
    limma_header   <- outputList[[i]][[3]]
    # create unique header for limma output
    colnames(limma_output) <- c(limma_lowLevel, paste0(limma_header[1, ], ".", limma_header[2, ]))
    
    # merge the outputs: data and limma output
    output <- dplyr::left_join(output, limma_output, by = limma_lowLevel)
    # merge the headers (headers and limma output) except the first column
    output_headers <- cbind(output_headers, limma_header)

  }
  
  return(list(headers = output_headers, data = output))
}


#
# MAIN
#

# Read YAML file
full_config <- read_yaml(opt$config)

# Merge LimmaCompare and General into a single list
config <- c(full_config$LimmaOneSample, full_config$General)

# override specific parameters from command-line arguments
config$report_infile <- opt$input
config$sample_table <- opt$samples
config$outfile <- opt$output

# prepare workspace
outdir <- dirname(config$outfile)  # Extract the directory path from the output file
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = TRUE)
}


config[['integrations']] <- c()
config[['lowLevel1']] <- c()
config[['lowLevel2']] <- c()

for (i in config[['ColumnNames']]) {
  config[['integrations']] <- c(config[['integrations']], i[[2]])
  config[['lowLevel1']] <- c(config[['lowLevel1']], i[[1]][1])
  config[['lowLevel2']] <- c(config[['lowLevel2']], i[[1]][2])
}

# Read with samples and parse it
sampleTable <- read.csv(config$sample_table, sep='\t', colClasses="character")
sampleGroups <- as.list(sampleTable)

for (i in names(sampleGroups)) {
  g <- c()
  for (j in sampleGroups[[i]]) {
    if (j!="") {
      g <- c(g, j)
    }
  }
  sampleGroups[[i]] <- g
}

# Set Logging file
logFile <- paste0(outdir, '/LimmaCompare.log')
basicConfig()
addHandler(writeToFile, logger='ReportLimma', file=logFile)


loginfo("Reading input files...", logger="ReportLimma")
headers <- read.csv(config$report_infile, header=FALSE, sep="\t", nrows=2)
data <- read.csv(config$report_infile, header=FALSE, sep="\t", skip=2)


loginfo("Calculating Limma pvalues...", logger="ReportLimma")
outputList <- calculatePvalues(headers, data, config, sampleGroups)


# Merge the list of outputs
loginfo("Merging limma results...", logger="ReportLimma")
results <- mergeOutputList(outputList, headers, data)


# Write output
loginfo("Writing results...", logger="ReportLimma")
# write the two header lines
write.table(results$headers, file = config$outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# append the data (without headers)
write.table(results$data, file = config$outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE, na = "")


loginfo("End script", logger="ReportLimma")
