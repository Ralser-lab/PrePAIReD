#' Function used to import DIA-NN into R
#' @param diann_report_file Path to the DIA-NN report file
#' @param norm_intensity Boolean; whether to use raw or normalized intensities from the DIA-NN report
#' @param fasta_file Path to the FASTA files used to annotate the protein IDs to genes
#' @param organism Name of organism; used to extract protein description from the FASTA file
#' @param metadata_file Path to the metadata file
#' @return List object with 3 data.tables: (1) `peptides` contains peptide intensities in all samples, (2) `proteins` contains protein IDs and gene names associated with the peptides, and (3) `samples` contains sample metadata
import_dataset <- function(diann_report_file, norm_intensity, fasta_file, organism = "HUMAN", metadata_file) {
  require(msdap)
  start.time <- Sys.time()
  # Import DIA-NN report file
  dataset <- import_dataset_diann(filename = diann_report_file, use_normalized_intensities = norm_intensity, do_plot = F)
  end.time <- Sys.time()
  print("Finished importing DIA-NN dataset")
  print(end.time - start.time)

  if (!is.null(fasta_file)) {
    # Annotate proteins in the data with gene names
    dataset <- import_fasta(dataset, files = fasta_file)
    if (!is.null(organism)) {
      # Extract gene descriptions from the FASTA entries
      dataset$proteins$description <- sapply(dataset$proteins$fasta_headers, function(header) {
        if (!grepl(";", header)) {
          substring(header, first = gregexpr(pattern = paste("_", toupper(organism), sep = ""), header)[[1]][1] + 7, last = gregexpr(pattern = "OS", header)[[1]][1] - 2)
        } else {
          paste(Reduce("c", lapply(strsplit(header, ";")[[1]], function(subheader) substring(subheader, first = gregexpr(pattern = paste("_", toupper(organism)), subheader)[[1]][1] + 7, last = gregexpr(pattern = "OS", subheader)[[1]][1] - 2))), collapse = ";")
        }
      })
    }
  }

  # Import sample metadata
  dataset <- import_sample_metadata(dataset, filename = metadata_file)
  return(dataset)
}

#' Function to filter out non-tryptic files
#' @param dataset List object with all DIA-NN data
#' @param tryptic_file File containing
tryptic_filt <- function(dataset, tryptic_file, venn_pdf) {
  tryptic_peps <- read.table(tryptic_file, header = T)
  tryptic <- unique(tryptic_peps$peptide[tryptic_peps$tryptic])

  # Determine how many proteins are tryptic and/or proteotypic
  proteotypic <- unique(dataset$peptides$sequence_plain[!grepl(";", dataset$peptides$protein_id)])

  pdf(venn_pdf)
  print(ggVennDiagram::ggVennDiagram(list(
    dataset = unique(dataset$peptides$sequence_plain),
    tryptic = tryptic[tryptic %in% unique(dataset$peptides$sequence_plain)],
    proteotypic = proteotypic
  )))
  graphics.off()

  # Remove non-tryptic peptides
  dataset <- remove_exclude_peps(dataset, unique(dataset$peptides$peptide_id[!(dataset$peptides$sequence_plain %in% tryptic)]))
  return(dataset)
}

# !/usr/bin/env Rscript
# if (!requireNamespace("optparse", quietly = TRUE)) {
#   install.packages("optparse")
# }
# if (!requireNamespace("yaml", quietly = TRUE)) {
#   install.packages("yaml")
# }

library(optparse)
library(yaml)
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
  make_option(c("-d", "--diann_file"),
    type = "character",
    help = "Path to the DIA-NN report file"
  ),
  make_option(c("-n", "--norm_intensity"),
    type = "logical", default = TRUE,
    help = "logical; whether to use raw or normalized intensities from the DIA-NN report [default= %default]"
  ),
  make_option(c("-f", "--fasta_file"),
    type = "character",
    help = "FASTA file used to annotate the proteins"
  ),
  make_option(c("-r", "--organism"),
    type = "character", default = "HUMAN",
    help = "Organism used to extract gene descriptions from FASTA file [default= %default]"
  ),
  make_option(c("-m", "--metadata_file"),
    type = "character",
    help = "Tab-delimited EXCEL file with metadata sample information"
  ),
  make_option(c("-t", "--tryptic_file"),
    type = "character",
    help = "Tab-separated text file with protein IDs and column `tryptic` on whether a protein is tryptic or not (TRUE/FALSE)"
  ),
  make_option(c("-o", "--output_dir"),
    type = "character", default = getwd(),
    help = "Path to save the import DIA-NN data into [default= %default]"
  ),
  make_option(c("-s", "--diann_stats"),
    type = "character", default = NULL,
    help = "Path to DIA-NN statistical report"
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
opt_parser <- OptionParser(option_list = option_list, description = "This R script is used to import DIA-NN data into R and save it as a RDS file with a `list` object.")
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
if (!exists("fasta_file")) fasta_file <- NULL
if (!exists("organism")) organism <- NULL
if (!exists("tryptic_file")) tryptic_file <- NULL
if (!exists("diann_stats")) diann_stats <- NULL

cmd.args <- commandArgs()
m <- regexpr("(?<=^--file=).+", cmd.args, perl = TRUE)
script.dir <- dirname(regmatches(cmd.args, m))
source(file.path(script.dir, "utils.R"))
cat("Parameters in use:\n")
cat("DIA-NN file =", file.path(input_dir, diann_file), "\n")
cat("Intensity used is normalized =", norm_intensity, "\n")
cat("FASTA file for protein annotatation =", fasta_file, "\n")
cat("Organism name to extract gene descriptions =", organism, "\n")
cat("Metadata file =", file.path(input_dir, metadata_file), "\n")
cat("Output directory =", file.path(input_dir, output_dir), "\n")

dir.create(file.path(input_dir, output_dir), showWarnings = F, recursive = T)
outFile <- file.path(input_dir, output_dir, "imported_data.rds")

# Import dataset
dataset <- import_dataset(
  diann_report_file = file.path(input_dir, diann_file),
  norm_intensity = norm_intensity, fasta_file = fasta_file, organism = organism, metadata_file = file.path(input_dir, metadata_file)
)

# Calculate dataset statistics (number of samples, peptides and proteins, as well as %age of missing values)
## Peptide matrix
mat <- dfToMat(dataset$peptides, peptide_id ~ sample_id, col = "intensity")
peptideNA <- 100 * sum(is.na(mat)) / length(mat)
## Protein matrix
mat <- dfToMat(dataset$peptides, protein_id ~ sample_id, col = "intensity")
proteinNA <- 100 * sum(mat == 0) / length(mat)
## Create table with statistics
overall_stats <- data.frame(orig = c(nrow(dataset$samples), length(unique(dataset$peptides$peptide_id)), length(unique(dataset$peptides$protein_id)), peptideNA, proteinNA))
rownames(overall_stats) <- c("nSamples", "nPeps", "nProts", "percPepNA", "percProtNA")

# If tryptic file is set, filter out non-tryptic peptides
if (!is.null(tryptic_file)) {
  cat("File with tryptic proteins =", tryptic_file, "\n")
  dataset <- tryptic_filt(dataset, tryptic_file, venn_pdf = file.path(input_dir, output_dir, "Venn_tryptic_proteotypic.pdf"))
}

# Calculate dataset statistics (number of samples, peptides and proteins, as well as %age of missing values)
## Peptide matrix
mat <- dfToMat(dataset$peptides, peptide_id ~ sample_id, col = "intensity")
peptideNA <- 100 * sum(is.na(mat)) / length(mat)
## Protein matrix
mat <- dfToMat(dataset$peptides, protein_id ~ sample_id, col = "intensity")
proteinNA <- 100 * sum(mat == 0) / length(mat)
## Create table with statistics
overall_stats$tryptic <- c(nrow(dataset$samples), length(unique(dataset$peptides$peptide_id)), length(unique(dataset$peptides$protein_id)), peptideNA, proteinNA)
overall_stats$tryptic_perc <- 100 * overall_stats$tryptic / overall_stats$orig

# If DIA-NN report file is set, generate QC report for it
if (!is.null(diann_stats)) {
  cat("Generating DIA-NN statistics report based on file:", file.path(input_dir, diann_stats), "\n")
  metadata <- dataset$samples[, c("file_name", batch, sample_type)]
  colnames(metadata) <- c("File.Name", "batch", "type")
  diann_stats_factors <- c("Total.Quantity", "MS1.Signal", "MS2.Signal", "FWHM.Scans", "FWHM.RT", "Median.Mass.Acc.MS1", "Median.Mass.Acc.MS1.Corrected", "Median.Mass.Acc.MS2", "Median.Mass.Acc.MS2.Corrected", "Normalisation.Instability", "Median.RT.Prediction.Acc", "Average.Peptide.Length", "Average.Peptide.Charge", "Average.Missed.Tryptic.Cleavages")
  dir.create(file.path(input_dir, output_dir, "QC_diann"), showWarnings = F, recursive = T)
  diann_stats <- qc_diann_figures(diann_stats_file = file.path(input_dir, diann_stats), metadata = metadata, factors = diann_stats_factors, outDir = file.path(input_dir, output_dir, "QC_diann"))
  write.xlsx(x = diann_stats, file = file.path(input_dir, output_dir, "QC_diann", "DIA-NN_stats.xlsx"))
}

# If stated as TRUE, generate QC report
if (qc_report) {
  cat("Generating a msdap QC report\n")
  cat("Sample type =", sample_type, "\n")
  cat("Batch =", batch, "\n")
  obj <- remove_exclude_samps(dataset, dataset$samples$sample_id[dataset$samples$exclude])
  dir.create(file.path(input_dir, output_dir, "QC_msdap"), showWarnings = F, recursive = T)
  msdap_stats <- qc_msdap_figures(dataset = obj, file.path(input_dir, output_dir, "QC_msdap"), intensity = "intensity", type = sample_type, batch = batch)
  write.xlsx(x = msdap_stats, file = file.path(input_dir, output_dir, "QC_msdap", "msdap_stats.xlsx"))
}

# Save dataset into file
saveRDS(object = dataset, file = outFile)
write.xlsx(x = overall_stats, file = file.path(input_dir, output_dir, "overall_stats.xlsx"))