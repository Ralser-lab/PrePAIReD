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