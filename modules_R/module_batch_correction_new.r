#!/usr/bin/env Rscript
#===========================================================
# Batch Correction Module (R)
# author: Roza Mizrak @rozsurmiz
# Compatible with existing helper_functions_R.R:
#   - requires: read_JSON(), run_with_logging()
#   - DOES NOT require: helper_init_module(), helper_close_module()
    #         ┌──────────────────────────────┐
    #         │    R Batch Correction Script │
    #         └──────────────────────────────┘
    #                      ▲
    #                      │
   #     ┌─────────────────┴─────────────────┐
  #      │                                   │
  #      │                                   │
 #Python wrapper                        Command-line call
 #(called from your                     (run Rscript directly)
 #cross-language pipeline)
    #    │                                   │
   #     ▼                                   ▼
  #┌─────────────────┐                 ┌─────────────┐
  #│ batchcorrect()  │◄────────────────┤  main()     │
  #│  (Public entry  │   function call  │  (CLI)      │
  #│   point)        │────────────────►│             │
  #└─────────────────┘   passes args    └─────────────┘
    #     │
   #      ▼
  #┌───────────────────────────┐
  #│ Core batch correction     │
  #│   logic:                  │
  #│   - read inputs           │
  #│   - align samples         │
  #│   - pick method (LIMMA,   │
  #│     ComBat, Harmony)      │
  #│   - write outputs         │
  #└───────────────────────────┘
#===========================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(limma)
  library(optparse)
  # Optional methods; loaded on demand
  suppressWarnings(suppressMessages(requireNamespace("sva", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("harmony", quietly = TRUE)))
})

# ---- Load your helpers (must define read_JSON(), run_with_logging()) ----
source("helper_functions_R.R")

# ---- small utils ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

.have_fun <- function(fname) exists(fname, mode = "function", inherits = TRUE)

.ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

.align_by_intersect <- function(data_df, meta_df) {
  # Expectation: rows = features, cols = samples
  common <- intersect(colnames(data_df), rownames(meta_df))
  drop_data <- setdiff(colnames(data_df), common)
  drop_meta <- setdiff(rownames(meta_df), common)

  warns <- character(0)
  if (length(drop_data)) {
    warns <- c(warns, paste0(
      "###WARNING### ", length(drop_data),
      " data column(s) not in metadata (dropped): ",
      paste(drop_data, collapse = ", ")
    ))
  }
  if (length(drop_meta)) {
    warns <- c(warns, paste0(
      "###WARNING### ", length(drop_meta),
      " metadata row(s) not in data (ignored): ",
      paste(drop_meta, collapse = ", ")
    ))
  }
  list(
    data_df = data_df[, common, drop = FALSE],
    meta_df = meta_df[common, , drop = FALSE],
    log     = if (length(warns)) paste(warns, collapse = "\n") else ""
  )
}

#===========================================================
# Core batch correction functions
#===========================================================

# limma::removeBatchEffect
batchcorrect_limma <- function(Data, batchER) {
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData

  batchF1 <- if (!is.null(batchER$batch))  factor(md[[batchER$batch]])  else NULL
  batchF2 <- if (!is.null(batchER$batch2)) factor(md[[batchER$batch2]]) else NULL
  covs    <- if (!is.null(batchER$covariates)) as.matrix(md[, batchER$covariates, drop = FALSE]) else NULL

  if (is.null(batchER$Class)) {
    corrected <- removeBatchEffect(DM, batch = batchF1, batch2 = batchF2, covariates = covs)
  } else {
    Class   <- factor(md[[batchER$Class]])
    designM <- model.matrix(~0 + Class)
    corrected <- removeBatchEffect(DM, batch = batchF1, batch2 = batchF2, covariates = covs, design = designM)
  }

  Data$Matrix <- corrected
  Data
}

# sva::ComBat (optional)
batchcorrect_combat <- function(Data, batchER) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required for method 'COMBAT' but is not installed.")
  }
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData

  par_prior <- batchER$par.prior %||% TRUE
  fitsamples <- batchER$fitsamples

  mat <- DM
  batch_vec <- md[[batchER$batch]]

  # Optional fit samples handling
  if (!is.null(fitsamples)) {
    fs <- intersect(fitsamples, rownames(md))
    if (length(fs) > 0) {
      # Ensure all batch levels present; if not, restrict to present ones
      present_levels <- unique(md[fs, batchER$batch])
      keep <- md[, batchER$batch] %in% present_levels
      if (!all(keep)) {
        message("###WARNING###: Restricting to batches present in fit samples for ComBat.")
        md  <- md[keep, , drop = FALSE]
        mat <- mat[, rownames(md), drop = FALSE]
        batch_vec <- md[[batchER$batch]]
      }
    }
  }

  mod <- NULL
  if (!is.null(batchER$covariates_preserve) && length(batchER$covariates_preserve)) {
    mod <- stats::model.matrix(~ ., data = md[, batchER$covariates_preserve, drop = FALSE])
  }

  corrected <- sva::ComBat(
    dat = as.matrix(mat),
    batch = as.factor(batch_vec),
    mod = mod,
    par.prior = par_prior,
    prior.plots = FALSE
  )

  # If we subset for fit, corrected already matches mat's columns
  Data$Matrix <- corrected
  Data
}

# harmony::HarmonyMatrix over PCA (optional)
batchcorrect_harmony <- function(Data, batchER) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' is required for method 'HARMONY' but is not installed.")
  }
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData
  n_pcs <- batchER$n_pcs %||% 20

  pca <- prcomp(t(DM), center = TRUE, scale. = FALSE)
  pcs <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]

  vars_use <- unique(na.omit(c(batchER$batch, batchER$covariates)))
  h_args   <- batchER$harmony_args %||% list()

  harmony_pcs <- do.call(
    harmony::HarmonyMatrix,
    c(list(data_mat = pcs, meta_data = md, vars_use = vars_use, do_pca = FALSE), h_args)
  )

  loadings <- pca$rotation[, seq_len(ncol(harmony_pcs)), drop = FALSE]
  centered <- t(DM) - matrix(pca$center, nrow = ncol(DM), ncol = nrow(DM), byrow = TRUE)
  corrected_centered <- harmony_pcs %*% t(loadings)
  corrected <- t(corrected_centered + matrix(pca$center, nrow = ncol(DM), ncol = nrow(DM), byrow = TRUE))

  Data$Matrix <- corrected
  Data
}

# Dispatcher
.dispatch_batch <- function(Data, batchER) {
  method <- toupper(batchER$method %||% "LIMMA")
  if (method == "LIMMA")   return(batchcorrect_limma(Data, batchER))
  if (method == "COMBAT")  return(batchcorrect_combat(Data, batchER))
  if (method == "HARMONY") return(batchcorrect_harmony(Data, batchER))
  stop("Unknown batch correction method: ", method)
}

#===========================================================
# Public entrypoint (wrapper calls this)
#===========================================================
batchcorrect <- function(data_df_fp,
                         meta_df_fp,
                         param_fp,
                         out_dir_fp,
                         mask_fp,
                         log_file_fp,
                         counter) {
  # ---- Local INIT (no helper_init_module assumed) ----
  if (!file.exists(data_df_fp) || !file.exists(meta_df_fp)) {
    stop("Input data or metadata file does not exist.")
  }
  data_df <- read.delim(data_df_fp, check.names = FALSE, row.names = 1)
  meta_df <- read.delim(meta_df_fp, check.names = FALSE, row.names = 1)
  params  <- read_JSON(param_fp)

  # Sample alignment
  ali <- .align_by_intersect(data_df, meta_df)
  data_df <- ali$data_df
  meta_df <- ali$meta_df
  log_string <- if (nzchar(ali$log)) paste0(ali$log, "\n") else ""

  # Params normalization (Python-like keys supported)
  method     <- toupper(params$method %||% "LIMMA")
  batch_cols <- params$batch_cols_meta %||% params$batch %||% "batch"
  if (length(batch_cols) < 1) stop("Missing batch column; expected 'batch' or 'batch_cols_meta' in params.")
  batch1 <- batch_cols[[1]]
  batch2 <- if (length(batch_cols) >= 2) batch_cols[[2]] else NULL

  batchER <- list(
    method      = method,
    batch       = batch1,
    batch2      = batch2,
    covariates  = params$covariates %||% params$combat_cols_adjust_meta %||% NULL,
    Class       = params$Class %||% NULL,
    par.prior   = params$par.prior %||% TRUE,
    n_pcs       = params$n_pcs %||% 20,
    harmony_args = params$harmony_args %||% list(),
    covariates_preserve = params$combat_cols_preserve_meta %||% NULL,
    fitsamples  = params$fitsamples_idx_list %||% NULL
  )

  Data <- list(Matrix = as.matrix(data_df), MetaData = meta_df)

  # ---- RUN ----
  corrected <- .dispatch_batch(Data, batchER)

  # ---- Local CLOSE (no helper_close_module assumed) ----
  .ensure_dir(out_dir_fp)
  suffix  <- sprintf("%s_%02d", format(Sys.time(), "%Y%m%d"), as.integer(counter %||% 1))
  outname <- sprintf("corrected_%s.tsv", suffix)
  out_fp  <- file.path(out_dir_fp, outname)
  write.table(corrected$Matrix, file = out_fp, sep = "\t", quote = FALSE, col.names = NA)

  if (!is.null(log_file_fp)) {
    cat(paste0("Batch correction complete. Output: ", out_fp, "\n"), file = log_file_fp, append = TRUE)
    if (nzchar(log_string)) cat(log_string, file = log_file_fp, append = TRUE)
  } else {
    message("Batch correction complete. Output written to: ", out_fp)
    if (nzchar(log_string)) message(log_string)
  }

  invisible(TRUE)
}

#===========================================================
# CLI
#===========================================================
main <- function() {
  option_list <- list(
    make_option(c("-d","--datapath"),  type = "character", help = "Data TSV filepath"),
    make_option(c("-m","--metapath"),  type = "character", help = "Metadata TSV filepath"),
    make_option(c("-p","--parampath"), type = "character", help = "JSON params filepath"),
    make_option(c("-o","--outpath"),   type = "character", help = "Output directory"),
    make_option(c("-k","--maskpath"),  type = "character", help = "Mask filepath (optional)"),
    make_option(c("-l","--log_file"),  type = "character", help = "Logfile filepath"),
    make_option(c("-c","--counter"),   type = "integer",   help = "Module counter")
  )
  opts <- parse_args(OptionParser(option_list = option_list))

  # Use your robust logger from helpers
  run_with_logging(
    expr = function() {
      batchcorrect(
        data_df_fp   = opts$datapath,
        meta_df_fp   = opts$metapath,
        param_fp     = opts$parampath,
        out_dir_fp   = opts$outpath,
        mask_fp      = opts$maskpath,
        log_file_fp  = opts$log_file,
        counter      = opts$counter
      )
    },
    log_path = opts$log_file
  )
}

if (sys.nframe() == 0) main()