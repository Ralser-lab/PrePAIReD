#!/usr/bin/env Rscript

##### Filtering samples #####
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

#####  #####

##### Filtering pepitdes #####
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
    help = "Path to directory with Rdata and overall statistics"
  ),
  make_option(c("-p", "--pep_thresh"),
    type = "numeric", default = 0.3,
    help = "Threshold to filter out low-presence peptides [default= %default]"
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
opt_parser <- OptionParser(option_list = option_list, description = "This R script is used to filter out low-presence peptides in the dataset `list` object.")
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

# Remove peptides expressed in too few samples
obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
mat <- dfToMat(obj$peptides, col = "intensity")
matFilt <- mat[rowSums(!is.na(mat)) >= ncol(mat) * pep_thresh, ]
dataset <- remove_exclude_peps(dataset, setdiff(rownames(mat), rownames(matFilt)))

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
mat <- dfToMat(dataset$peptides, peptide_id ~ sample_id, col = "intensity_by_group")
nSamps <- ncol(mat)
nPeps <- nrow(mat)
peptideNA <- 100 * sum(is.na(mat)) / length(mat)
## Protein matrix
mat <- dfToMat(dataset$peptides, protein_id ~ sample_id, col = "intensity_by_group")
proteinNA <- 100 * sum(mat == 0) / length(mat)
nProts <- nrow(mat)
## Create table with statistics
overall_stats <- read.xlsx(file.path(input_dir, data_dir, "overall_stats.xlsx"))
overall_stats$filtPeps <- c(nSamps, nPeps, nProts, peptideNA, proteinNA)
overall_stats$filtPeps_perc <- 100 * overall_stats$filtPeps / overall_stats$orig

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
outFile <- file.path(input_dir, output_dir, "filtPeps_data.rds")
dir.create(file.path(input_dir, output_dir), showWarnings = F, recursive = T)
saveRDS(object = dataset, file = outFile)
write.xlsx(x = overall_stats, file = file.path(input_dir, output_dir, "overall_stats.xlsx"))

##### #####

##### Imputation #####
# Support function to retrieve level of missingness in rows
get_missingness <- function(data) {
  rowMeans(is.na(data)) * 100 # Get percentage of missing values per peptide (rows)
}

# !/usr/bin/env Rscript
library(optparse)
library(yaml)
library(msdap)
library(openxlsx)
library(reshape2)

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
    help = "Path to directory with Rdata file and overall statistics"
  ),
  make_option(c("-f", "--impute_factor"),
    type = "character", default=NULL,
    help = "Factor used to impute data"
  ),
  make_option(c("-m", "--impute_method"),
    type = "character", default = "knn",
    help = "Method used to impute the data [default= %default]"
  ),
  make_option(c("-t", "--miss_thresh"),
    type = "numeric",
    help = "Missingness threshold to filter peptides to be imputed"
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
opt_parser <- OptionParser(option_list = option_list, description = "This R script is used to impute the dataset `list` object.")
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

# Create intensity matrix
obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
mat <- dfToMat(obj$peptides, peptide_id ~ sample_id, col = "intensity_by_group")
data_imputed <- mat

if (!is.null(impute_factor)) {
  # Impute data for specific factor
  impute_groups <- unique(obj$samples[[impute_factor]])

  for (impute_g in impute_groups) {
    subset_samples <- obj$samples$sample_id[obj$samples[[impute_factor]] == impute_g]
    subset_data <- mat[, colnames(mat) %in% subset_samples] # Extract data for the factor to impute

    missingness <- get_missingness(subset_data)
    peptides_to_impute <- names(missingness[missingness <= miss_thresh])

    if (length(peptides_to_impute) > 0) {
      start.time <- Sys.time()
      tmp <- impute.matrix(subset_data[peptides_to_impute, ], method = impute_method)
      end.time <- Sys.time()
      cat("Finished imputing element", impute_g, "in imputation factor", impute_factor,"\n")
      print(end.time - start.time)
      data_imputed[match(rownames(tmp), rownames(data_imputed)), match(colnames(tmp), colnames(data_imputed))] <- tmp
    }
  }
} else {
  missingness <- get_missingness(mat)
  peptides_to_impute <- names(missingness[missingness <= miss_thresh])

  if (length(peptides_to_impute) > 0) {
    start.time <- Sys.time()
    tmp <- impute.matrix(mat[peptides_to_impute, ], method = impute_method)
    end.time <- Sys.time()
    cat("Finished imputing factor", impute_g, "\n")
    print(end.time - start.time)
    data_imputed[match(rownames(tmp), rownames(data_imputed)), match(colnames(tmp), colnames(data_imputed))] <- tmp
  }
}
df <- melt(data_imputed)
colnames(df) <- c("peptide_id", "sample_id", "intensity_by_group_imputed")
dataset$peptides <- merge(dataset$peptides, df, all = T)

# Imputation resulted in some rows with NAs; merge the peptide and sample information back into the data 
dataset$peptides$isdecoy <- FALSE
dataset$peptides$detect <- TRUE
## Peptide info
tmp <- unique(dataset$peptides[!is.na(dataset$peptides$sequence_modified), c("peptide_id", "protein_id", "sequence_plain", "sequence_modified")])
dataset$peptides$protein_id <- dataset$peptides$sequence_modified <- dataset$peptides$sequence_plain <- NULL
dataset$peptides <- merge(tmp, dataset$peptides)
## Sample info
tmp <- unique(dataset$peptides[!is.na(dataset$peptides$shortname), c("shortname", "sample_id")])
dataset$peptides$shortname <- NULL
dataset$peptides <- merge(tmp, dataset$peptides)

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
overall_stats <- read.xlsx(file.path(input_dir, overall_stats))
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
outFile <- file.path(input_dir, output_dir, "imputed_data.rds")
dir.create(file.path(input_dir, output_dir), showWarnings = F, recursive = T)
saveRDS(object = dataset, file = outFile)
write.xlsx(x = overall_stats, file = file.path(input_dir, output_dir, "overall_stats.xlsx"))

##### #####

##### Batch correction #####
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
    help = "Path to directory with Rdata file and overall statistics"
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
opt_parser <- OptionParser(option_list = option_list, description = "This R script is used to impute the dataset `list` object.")
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

obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
obj$peptides <- remove_batch_effect(obj, batch = batch, col = "intensity_by_group")
outDir <- file.path(input_dir, output_dir)
dir.create(outDir, showWarnings = F, recursive = T)

qc_df_batchCorr <- qc_msdap_figures(obj, outDir, batch = batch, type = sample_type, intensity = "intensity_by_group")

cat("Calculating PCA for plate-corrected data\n")
obj$samples$ms_batch_n <- as.numeric(gsub("b", "", obj$samples$ms_batch))
obj$samples$analysis_plate96_n <- as.numeric(gsub("P", "", obj$samples$analysis_plate96))

mat <- dfToMat(obj$peptides, peptide_id ~ sample_id, col = "intensity_by_group")

res <- pca_msdap_res(dataset = obj, impute = "ppca", features = "peptides", intensity = "intensity_by_group", quanti_vars = c("detected_peptides", "detected_proteins", "ms_batch_n", "analysis_plate96_n"), quali_vars = c("analysis_plate96", "ms_batch"), calc.nPCs = 15)
plot.nPCs <- pca_optimal(res$res.pca)

dir.create(file.path(outDir, "PCA"), showWarnings = F, recursive = T)
pca_msdap_figures(res.pca = res$res.pca, outDir = file.path(outDir, "PCA"), type = "analysis_plate96", calc.nPCs = res$calc.nPCs, plot.nPCs = plot.nPCs)
##### #####
