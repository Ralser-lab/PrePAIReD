#!/usr/bin/env Rscript

library(optparse)
library(yaml)
library(msdap)
library(openxlsx)

option_list <- list(
  make_option(c("-c", "--config"),
    type = "character",
    help = "Path to configuration file (YAML format)"
  ),
  make_option(c("-i", "--input_dir"),
    type = "character",
    help = "Path to main input directory"
  ),
  make_option(c("-d", "--data_dir"),
    type = "character",
    help = "Path to Rdata file and overall statistics"
  ),
  make_option(c("-t", "--thresh"),
    type = "numeric",
    help = "Threshold to filter out outlier samples"
  ),
  make_option(c("-n", "--norm_methods"),
    type = "character", default = c("loess", "modebetween_protein"),
    help = "Normalization methods used for the data [default= %default]"
  ),
  make_option(c("-r", "--rollup_algorithm"),
    type = "character", default = "maxlfq_diann",
    help = "Algorithm used to rollup peptide to protein levels [default= %default]"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character", default = getwd(),
    help = "Path to save the import DIA-NN data into [default= %default]"
  ),
  make_option(c("-q", "--qc_report"),
    type = "logical", default = TRUE,
    help = "Whether to generate the msdap QC report or not [default= %default]"
  ),
  make_option(c("-b", "--batch"),
    type = "character", default = "ms_batch",
    help = "If QC report is generated, set batches to be plotted [default= %default]"
  ),
  make_option(c("-y", "--sample_type"),
    type = "character", default = "sample_type",
    help = "If QC report is generated, set sample type to be plotted [default= %default]"
  )
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list, description = "This R script is used to filter out outlier samples with too few peptides measured in the dataset `list` object.")
opt <- parse_args(opt_parser)

# Function to load config file
load_config <- function(config_path) {
  config <- yaml::read_yaml(config_path)
  return(config)
}

# Load parameters from config or command line
if (!is.null(opt$config)) {
  config <- load_config(opt$config)
  for (i in 1:length(config)) assign(x = names(config)[i], value = config[[i]])
  for (i in 1:length(opt)) if (!exists(names(opt)[i])) assign(x = names(opt)[i], value = opt[[i]])
} else {
  for (i in 1:length(opt)) assign(x = names(opt)[i], value = opt[[i]])
}

cmd.args <- commandArgs()
m <- regexpr("(?<=^--file=).+", cmd.args, perl = TRUE)
script.dir <- dirname(regmatches(cmd.args, m))
source(file.path(script.dir, "utils.R"))
# cat("Parameters in use:\n")
# cat("DIA-NN file =", file.path(input_dir, diann_file), "\n")
# cat("Intensity used is normalized =", norm_intensity, "\n")
# cat("FASTA file for protein annotatation =", fasta_file, "\n")
# cat("Organism name to extract gene descriptions =", organism, "\n")
# cat("Metadata file =", file.path(input_dir, metadata_file), "\n")
# cat("Output directory =", file.path(input_dir, output_dir), "\n")

# Read data file from previous pre-processing step
dataset <- readRDS(list.files(file.path(input_dir, data_dir), pattern = "rds", full.names = T))

# Remove outlier samples with too few peptides measured
mat <- dfToMat(dataset$peptides, col = "intensity")
dataset$samples$percNA <- sapply(dataset$samples$sample_id, function(s) {
  tmp <- mat[, colnames(mat) %in% s]
  sum(is.na(tmp)) / length(tmp)
})
dataset$samples$exclude <- ifelse(dataset$samples$percNA < thresh, FALSE, TRUE)

# Remove dt_pep_group and groups sections in dataset, and set them again
dataset$dt_pep_group <- dataset$groups <- NULL
dataset <- cache_filtering_data_noCoV(dataset)
# Normalize data
cat("Normalising data\n")
dataset <- filter_dataset(
  dataset = dataset,
  filter_min_detect = 3, # each peptide must have a good confidence score in at least N samples per group
  filter_min_quant = 3, # similarly, the number of reps where the peptide must have a quantitative value
  by_group = T, all_group = F, by_contrast = F,
  rollup_algorithm = rollup_algorithm, norm_algorithm = norm_methods
) # loess calls the cyclic loess function normalizeCyclicLoess from limma (recommended by Vadim), and the modebetween_protein scales between all pairwise sample groups on a protein level, such that the log-fold-change mode is zero

# Calculate dataset statistics (number of samples, peptides and proteins, as well as %age of missing values)
## Peptide matrix
obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
mat <- dfToMat(obj$peptides, peptide_id ~ sample_id, col = "intensity_by_group")
nSamps <- ncol(mat)
nPeps <- nrow(mat)
peptideNA <- 100 * sum(is.na(mat)) / length(mat)
## Protein matrix
mat <- dfToMat(obj$peptides, protein_id ~ sample_id, col = "intensity_by_group")
proteinNA <- 100 * sum(mat == 0) / length(mat)
nProts <- nrow(mat)
## Create table with statistics
overall_stats <- read.xlsx(file.path(input_dir, data_dir, "overall_stats.xlsx"))
overall_stats$filtSamps <- c(nSamps, nPeps, nProts, peptideNA, proteinNA)
overall_stats$filtSamps_perc <- 100 * overall_stats$filtSamps / overall_stats$orig

# If stated as TRUE, generate QC report
if (qc_report) {
  cat("Generating a msdap QC report\n")
  cat("Sample type =", sample_type, "\n")
  cat("Batch =", batch, "\n")
  obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
  dir.create(file.path(input_dir, output_dir, "QC_msdap"), showWarnings = F, recursive = T)
  msdap_stats <- qc_msdap_figures(dataset = obj, file.path(input_dir, output_dir, "QC_msdap"), intensity = "intensity_by_group", type = sample_type, batch = batch)
  write.xlsx(x = msdap_stats, file = file.path(input_dir, output_dir, "QC_msdap", "msdap_stats.xlsx"))
}

# Create directory for QC results
outFile <- file.path(input_dir, output_dir, "filtSamps_data.rds")
dir.create(file.path(input_dir, output_dir), showWarnings = F, recursive = T)
saveRDS(object = dataset, file = outFile)
write.xlsx(x = overall_stats, file = file.path(input_dir, output_dir, "overall_stats.xlsx"))