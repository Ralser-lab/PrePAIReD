#===========================================================
# Batch Correction Module (R)
# author: Roza Mizrak @rozsurmiz
#===========================================================

#===========================================================
# Core batch correction functions
#===========================================================

#' Batch correction via limma::removeBatchEffect
batchcorrect_limma <- function(Data, batchER) {
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData
  batchF1 <- if (!is.null(batchER$batch)) factor(md[[batchER$batch]]) else NULL
  batchF2 <- if (!is.null(batchER$batch2)) factor(md[[batchER$batch2]]) else NULL
  covs    <- if (!is.null(batchER$covariates)) as.matrix(md[, batchER$covariates, drop=FALSE]) else NULL
  if (is.null(batchER$Class)) {
    corrected <- removeBatchEffect(DM, batch=batchF1, batch2=batchF2, covariates=covs)
  } else {
    Class   <- factor(md[[batchER$Class]]);
    designM <- model.matrix(~0 + Class)
    corrected <- removeBatchEffect(DM, batch=batchF1, batch2=batchF2,
                                   covariates=covs, design=designM)
  }
  Data$Matrix <- corrected; return(Data)
}

#' Batch correction via sva::ComBat
batchcorrect_combat <- function(Data, batchER) {
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData
  par_prior <- if (!is.null(batchER$par.prior)) batchER$par.prior else TRUE
  corrected <- ComBat(dat=DM, batch=md[[batchER$batch]], mod=NULL,
                      par.prior=par_prior, prior.plots=FALSE)
  Data$Matrix <- corrected; return(Data)
}

#' Batch correction via harmony integration
batchcorrect_harmony <- function(Data, batchER) {
  set.seed(42)
  DM <- Data$Matrix; md <- Data$MetaData
  n_pcs      <- if (!is.null(batchER$n_pcs)) batchER$n_pcs else 20
  pca <- prcomp(t(DM), center=TRUE, scale.=FALSE)
  pcs <- pca$x[, seq_len(n_pcs), drop=FALSE]
  extras <- batchER[setdiff(names(batchER), c("method","batch","batch2","covariates","Class","n_pcs","par.prior"))]
  harmony_pcs <- do.call(HarmonyMatrix, c(
    list(data_mat=pcs,
         meta_data=md,
         vars_use=c(batchER$batch, batchER$covariates),
         do_pca=FALSE),
    extras
  ))
  loadings <- pca$rotation[, seq_len(n_pcs), drop=FALSE]
  centered <- t(DM) - matrix(pca$center,
                             nrow=ncol(DM),
                             ncol=nrow(DM),
                             byrow=TRUE)
  corrected_centered <- harmony_pcs %*% t(loadings)
  corrected <- t(corrected_centered +
                 matrix(pca$center,
                        nrow=ncol(DM),
                        ncol=nrow(DM),
                        byrow=TRUE))
  Data$Matrix <- corrected; return(Data)
}

#' Dispatcher for batch correction methods
#batchcorrect <- function(Data, batchER) {
#  method <- toupper(batchER$method)
#  if (method == "LIMMA") return(batchcorrect_limma(Data, batchER))
#  if (method == "COMBAT") return(batchcorrect_combat(Data, batchER))
#  if (method == "HARMONY") return(batchcorrect_harmony(Data, batchER))
#  stop("Unknown method: ", batchER$method)
#}

#' Run batch correction from R (callable function)
#'
#' @param data_df_fp  filepath to the data matrix TSV
#' @param meta_df_fp  filepath to the metadata TSV
#' @param param_fp    filepath to the JSON parameters
#' @param out_dir_fp  directory for output files
#' @param mask_fp     filepath to mask file (currently unused)
#' @param log_file_fp filepath to write logs
#' @param counter     integer module counter
#' @export
run_batch_correction <- function(data_df_fp, meta_df_fp, param_fp,
                                 out_dir_fp, mask_fp, log_file_fp,
                                 counter) {
  # Setup logging environment
  assign("log_file_fp", log_file_fp, envir = environment(log_info))
  log_info(sprintf("### MODULE %s START ###", counter))

  # Load and validate input files
  if (!file.exists(data_df_fp) || !file.exists(meta_df_fp)) {
    log_error("Input data or metadata file does not exist")
  }
  data_df <- tryCatch(
    read.csv(data_df_fp, sep="	", row.names=1, check.names=FALSE),
    error = function(e) log_error(paste("Failed to read data_df_fp:", e$message))
  )
  meta_df <- tryCatch(
    read.csv(meta_df_fp, sep="	", row.names=1, check.names=FALSE),
    error = function(e) log_error(paste("Failed to read meta_df_fp:", e$message))
  )

  # Handle missing batch annotations
  batch_col <- "batch"
  if (!batch_col %in% colnames(meta_df)) {
    log_error(sprintf("Missing '%s' column in meta_df", batch_col))
  }
  missing_batch <- rownames(meta_df)[is.na(meta_df[[batch_col]]) | meta_df[[batch_col]] == ""]
  if (length(missing_batch) > 0) {
    log_warning(paste("Samples with missing batch:", paste(missing_batch, collapse=", ")))
    data_df <- data_df[, !colnames(data_df) %in% missing_batch, drop=FALSE]
    meta_df <- meta_df[!rownames(meta_df) %in% missing_batch, , drop=FALSE]
  }

  # Read JSON parameters and dynamic dispatch
  param_dict <- tryCatch(
    jsonlite::fromJSON(param_fp),
    error = function(e) log_error(paste("Failed to parse param_fp:", e$message))
  )
  # Determine method and function
  method   <- toupper(param_dict$method %||% "LIMMA")
  fn_name  <- paste0("batchcorrect_", tolower(method))
  if (!exists(fn_name, mode = "function")) {
    log_error(paste("Unknown batch correction method:", method))
  }
  fn        <- get(fn_name)
  fn_params <- names(formals(fn))
  # Build batchER list with only relevant params
  batchER  <- c(
    list(method = method, batch = batch_col),
    param_dict[names(param_dict) %in% fn_params]
  )

  # Prepare Data list
  data_arr  <- as.matrix(data_df)
  Data_list <- list(Matrix = data_arr, MetaData = meta_df)

  # Invoke batch correction
  corrected_Data_list <- tryCatch(
    do.call(fn, list(Data = Data_list, batchER = batchER)),
    error = function(e) log_error(paste("Batch correction failed:", e$message))
  )

  # Write output
  if (!dir.exists(out_dir_fp)) dir.create(out_dir_fp, recursive = TRUE)
  output_fp <- file.path(out_dir_fp, sprintf("corrected_%s.tsv", counter))
  tryCatch(
    write.table(corrected_Data_list$Matrix,
                file = output_fp,
                sep = "	",
                quote = FALSE,
                col.names = NA),
    error = function(e) log_error(paste("Failed to write output_fp:", e$message))
  )

  log_info(sprintf("### MODULE %s COMPLETED ###", counter))
}

#' Main entry for command-line execution
main <- function() {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option(c("-d","--datapath"), type="character", help="Data TSV filepath"),
    make_option(c("-md","--metapath"), type="character", help="Metadata TSV filepath"),
    make_option(c("-p","--parampath"), type="character", help="JSON params filepath"),
    make_option(c("-o","--outpath"), type="character", help="Output directory"),
    make_option(c("-ma","--maskpath"), type="character", help="Mask filepath (unused)"),
    make_option(c("-l","--log_file"), type="character", help="Logfile filepath"),
    make_option(c("-c","--counter"), type="integer",   help="Module counter")
  )
  opts <- parse_args(OptionParser(option_list = option_list))
  run_batch_correction(
    data_df_fp  = opts$datapath,
    meta_df_fp  = opts$metapath,
    param_fp    = opts$parampath,
    out_dir_fp  = opts$outpath,
    mask_fp     = opts$maskpath,
    log_file_fp = opts$log_file,
    counter     = opts$counter
  )
}

# Invoke when run as script
if (sys.nframe() == 0) main()
