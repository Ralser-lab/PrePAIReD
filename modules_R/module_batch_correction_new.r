#===========================================================
# Batch Correction Module (R)
# author: Roza Mizrak @rozsurmiz
#===========================================================
# Load helpers
library(jsonlite)
library(limma)

# now source the helpers by absolute path
source("helper_functions_R.R")
# defines read_JSON() and run_with_logging()

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
#' @export
#===========================================================
# The workhorse: read inputs, dispatch, write output
#===========================================================
batchcorrect <- function(data_df_fp,
                          meta_df_fp,
                          param_fp,
                          out_dir_fp,
                          mask_fp,
                          log_file_fp,
                          counter) {

  # generate run suffix from timestamp / counter
  suffix <- sprintf("%s_%02d", format(Sys.time(), "%Y%m%d"), counter)

  # 1) Load & validate input files
  if (!file.exists(data_df_fp) || !file.exists(meta_df_fp)) {
    stop("Input data or metadata file does not exist")
  }
  data_df <- read.delim(data_df_fp, check.names=FALSE, row.names=1)
  meta_df <- read.delim(meta_df_fp, check.names=FALSE, row.names=1)

  # 2) Drop samples missing batch
  if (!"batch" %in% colnames(meta_df)) {
    stop("Missing 'batch' column in metadata")
  }
  nas <- is.na(meta_df$batch) | meta_df$batch==""
  if (any(nas)) {
    warning("Dropping samples with missing batch: ",
            paste(rownames(meta_df)[nas], collapse=", "))
    data_df <- data_df[, !colnames(data_df) %in% rownames(meta_df)[nas], drop=FALSE]
    meta_df <- meta_df[!nas, , drop=FALSE]
  }

  # 3) Read params via read_JSON()
  params <- read_JSON(param_fp)
  method <- toupper(params$method %||% "LIMMA")
  fn_name <- paste0("batchcorrect_", tolower(method))
  if (!exists(fn_name, mode="function")) {
    stop("Unknown batch correction method: ", method)
  }
  fn <- get(fn_name)

  # Build batchER list
  valid_args <- names(formals(fn))
  batchER <- c(
    list(method=method, batch="batch"),
    params[names(params) %in% valid_args]
  )

  # 4) Run correction
  Data_list <- list(Matrix=as.matrix(data_df), MetaData=meta_df)
  corrected <- do.call(fn, list(Data=Data_list, batchER=batchER))

  # 5) Write output
  if (!dir.exists(out_dir_fp)) dir.create(out_dir_fp, recursive=TRUE)
  out_fp <- file.path(out_dir_fp, sprintf("corrected_%s.tsv", suffix))
  write.table(corrected$Matrix, file=out_fp,
              sep="\t", quote=FALSE, col.names=NA)

  message("Batch correction complete. Output written to: ", out_fp)
}

#===========================================================
# Main: parse CLI args and wrap with run_with_logging()
#===========================================================
main <- function() {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option(c("-d","--datapath"),   type="character", help="Data TSV filepath"),
    make_option(c("-md","--metapath"),   type="character", help="Metadata TSV filepath"),
    make_option(c("-p","--parampath"),   type="character", help="JSON params filepath"),
    make_option(c("-o","--outpath"),     type="character", help="Output directory"),
    make_option(c("-ma","--maskpath"),   type="character", help="Mask filepath (unused)"),
    make_option(c("-l","--log_file"),    type="character", help="Logfile filepath"),
    make_option(c("-c","--counter"),     type="integer",   help="Module counter")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  # Wrap the entire pipeline in your robust logger:
  run_with_logging(
    expr = function() {
      process_batch(
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

# If invoked directly, fire off main()
if (sys.nframe() == 0) main()