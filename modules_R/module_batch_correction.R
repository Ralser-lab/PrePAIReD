#!/usr/bin/env Rscript
#===========================================================
# Batch Correction Module (R) — refactored to follow the new CLI pattern
# author: Roza Mizrak @rozsurmiz 
# Uses helper_functions_R.R: make_parser(), validate_opts(), run_with_logging(), read_JSON()
# - Pure function: batch_correct(data, meta, params) -> corrected data.frame
# - CLI main(): parses args in Python-like style and handles I/O + logging

#         ┌──────────────────────────────┐
#         │    R Batch Correction Script │
#         └──────────────────────────────┘
#                      ▲
#                      │
#     ┌─────────────────┴─────────────────┐
#      │                                   │
#      │                                   │
#Python wrapper                        Command-line call
#(called from                        (run Rscript directly)
#cross-language pipeline)
#    │                                   │
#     ▼                                   ▼
#┌─────────────────┐                 ┌─────────────┐
#│ batchcorrect()  │◄────────────────┤  main()     │
#│  (Public entry  │   function call │  (CLI)      │
#│   point)        │────────────────►│             │
#└─────────────────┘   passes args   └─────────────┘
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
  library(optparse)
  library(data.table)
  suppressWarnings({ 
    # load commonly used methods if present
    if (requireNamespace("sva", quietly = TRUE)) library(sva)
    if (requireNamespace("limma", quietly = TRUE)) library(limma)
  })
})

suppressPackageStartupMessages({
  library(jsonlite)
  library(limma)
  library(optparse)
  # Optional methods; loaded on demand
  suppressWarnings(suppressMessages(requireNamespace("sva", quietly = TRUE)))
  suppressWarnings(suppressMessages(requireNamespace("harmony", quietly = TRUE)))
})

# ---- Load your helpers (must define read_JSON(), run_with_logging()) ----
# ---- Portable helper loader (no hard-coded paths) --------------------------
.load_helpers <- function() {
  env_path <- Sys.getenv("PIPELINE_HELPERS_R", unset = "")
  if (nzchar(env_path) && file.exists(env_path)) {
    source(env_path, chdir = TRUE); return(invisible(TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- tryCatch({
    f <- args[grep(file_arg, args, fixed = TRUE)]
    if (length(f)) normalizePath(sub(file_arg, "", f[1])) else NA_character_
  }, error = function(e) NA_character_)
  candidates <- unique(na.omit(c(
    if (!is.na(script_path)) file.path(dirname(script_path), "helper_functions_R.R") else NA_character_,
    if (!is.na(script_path)) file.path(dirname(script_path), "..", "helper_functions_R.R") else NA_character_,
    file.path(getwd(), "helper_functions_R.R")
  )))
  for (cand in candidates) {
    if (file.exists(cand)) { source(cand, chdir = TRUE); return(invisible(TRUE)) }
  }
  stop(
    paste0(
      "helper_functions_R.R not found.\n",
      "Set env var PIPELINE_HELPERS_R=/absolute/path/to/helper_functions_R.R\n",
      "or place helper_functions_R.R next to this script.\n",
      "Searched:\n - ", paste(candidates, collapse = "\n - ")
    ),
    call. = FALSE
  )
}
.load_helpers()
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

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
    stop("Package 'sva' is required for method 'COMBAT_R' but is not installed.")
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
    stop("Package 'harmony' is required for method 'HARMONY_R' but is not installed.")
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
  method <- toupper(batchER$correct_/Users/rozsurmiz/Documents/GitHub/Preprocessing/modules_R/test_batch_correction.shmethod %||% "LIMMA__R")
  if (method == "LIMMA_R")   return(batchcorrect_limma(Data, batchER))
  if (method == "COMBAT_R")  return(batchcorrect_combat(Data, batchER))
  if (method == "HARMONY_R") return(batchcorrect_harmony(Data, batchER))
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
  method     <- toupper(params$correct_method %||% "LIMMA_R")
  batch_cols <- params$batch_cols_meta %||% params$batch %||% "batch"
  if (length(batch_cols) < 1) stop("Missing batch column; expected 'batch' or 'batch_cols_meta' in params.")
  batch1 <- batch_cols[[1]]
  batch2 <- if (length(batch_cols) >= 2) batch_cols[[2]] else NULL

  batchER <- list(
    correct_method = method,
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
# Pure entrypoint (no I/O): returns corrected data.frame
#===========================================================
batch_correct <- function(data_df, meta_df, params) {
  # Align samples
  ali <- .align_by_intersect(data_df, meta_df)
  data_df <- ali$data_df
  meta_df <- ali$meta_df
  log_string <- if (nzchar(ali$log)) paste0(ali$log, "\n") else ""

  method     <- toupper(params$correct_method %||% "LIMMA_R")
  batch_cols <- params$batch_cols_meta %||% params$batch %||% "batch"
  if (length(batch_cols) < 1) stop("Missing batch column; expected 'batch' or 'batch_cols_meta' in params.")
  batch1 <- batch_cols[[1]]
  batch2 <- if (length(batch_cols) >= 2) batch_cols[[2]] else NULL

  batchER <- list(
    correct_method = method,
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
  corrected <- .dispatch_batch(Data, batchER)
  M_out <- corrected$Matrix
  # Return as data.frame preserving rownames
  out_df <- as.data.frame(M_out, check.names = FALSE)
  return(out_df)
}


#===========================================================
# CLI
#===========================================================
main <- function() {
  # Parse args using the shared helper
  parser <- make_parser()
  args   <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))
  validate_opts(args)

  # Read parameters
  params <- read_JSON(args$params)

  # Build input paths
  input_file <- file.path(args$data, sprintf("data_%02d.tsv", as.integer(args$counter) - 1))
  if (!file.exists(input_file)) {
    stop(sprintf("Input file not found: %s", input_file), call. = FALSE)
  }
  if (is.null(args$metadata) || !file.exists(args$metadata)) {
    stop("Batch correction requires a metadata file (--metadata).", call. = FALSE)
  }

  # ---- Load matrices with proper rownames ----
  DM_df <- data.table::fread(input_file, sep = "\t", data.table = FALSE, check.names = FALSE)
  rownames(DM_df) <- DM_df[[1]]; DM_df[[1]] <- NULL
  DM_df <- as.data.frame(DM_df, check.names = FALSE)
  
  meta_df <- data.table::fread(args$metadata, sep = "\t", data.table = FALSE, check.names = FALSE)
  rownames(meta_df) <- meta_df[[1]]; meta_df[[1]] <- NULL
  meta_df <- as.data.frame(meta_df, check.names = FALSE)
  
  # ---- Helpers for name reconciliation ----
  .normalize_runs <- function(x) {
    x <- gsub("\\\\", "/", x)              # Windows -> POSIX
    x <- sub(".*/", "", x)                 # basename
    # strip common vendor suffixes/extensions
    x <- sub("\\.rawIdx\\.wiff\\.dia$", "", x, ignore.case = TRUE)
    x <- sub("\\.wiff\\.dia$", "", x, ignore.case = TRUE)
    x <- sub("\\.raw$", "", x, ignore.case = TRUE)
    x <- sub("\\.mzML$", "", x, ignore.case = TRUE)
    x <- sub("\\.d$", "", x, ignore.case = TRUE)
    x
  }
  choose_meta_id <- function(meta, run_names_norm) {
    best <- NULL; best_score <- -1
    for (cn in names(meta)) {
      cand <- .normalize_runs(meta[[cn]])
      sc <- length(intersect(run_names_norm, cand))
      if (sc > best_score) { best <- cn; best_score <- sc }
    }
    list(best = best, score = best_score)
  }
  
  # ---- Normalize current data colnames up front ----
  colnames(DM_df) <- .normalize_runs(colnames(DM_df))
  
  # ---- Apply --ids mapping (normalize BOTH keys and values) ----
  if (!is.null(args$ids) && file.exists(args$ids)) {
    idmap <- data.table::fread(args$ids, sep = "\t", data.table = FALSE, check.names = FALSE)
    lnms <- tolower(names(idmap))
    # try by header names
    from_idx <- which(lnms %in% c("raw","raw_file","file","path","run","original","old","source"))
    to_idx   <- which(lnms %in% c("sample","sample.name","name","new","target","id"))
    # fallback by content (paths/extensions in 'from'; overlap to runs in 'to')
    if (length(from_idx) == 0)
      from_idx <- which(sapply(idmap, function(col) any(grepl("[/\\\\]|\\.raw|\\.wiff|\\.dia|\\.mzML", col, ignore.case=TRUE))))
    if (length(to_idx) == 0)
      to_idx <- which(sapply(idmap, function(col) length(intersect(.normalize_runs(col), colnames(DM_df))) > 0))
    if (length(from_idx) >= 1 && length(to_idx) >= 1) {
      from <- names(idmap)[from_idx[1]]; to <- names(idmap)[to_idx[1]]
      m <- setNames(.normalize_runs(idmap[[to]]), .normalize_runs(idmap[[from]]))
      # rename any columns that match a mapping key
      current <- colnames(DM_df)
      hit <- intersect(current, names(m))
      if (length(hit)) {
        idx <- match(hit, current)
        new <- current
        new[idx] <- unname(m[hit])
        colnames(DM_df) <- new
        message("###DEBUG### applied ids mapping: ", length(hit), " columns renamed via --ids (", from, " -> ", to, ")")
      } else {
        message("###WARNING### --ids provided but no keys matched current column names after normalization.")
      }
    } else {
      message("###WARNING### Could not auto-detect from/to columns in --ids file; skipping ids mapping.")
    }
  }
  
  
  # 5) **Normalize metadata rownames**  <-- ADD THIS LINE
  rownames(meta_df) <- .normalize_runs(rownames(meta_df))
  
  # ---- If metadata rownames aren’t sample IDs, pick the best column and use it ----
  run_names <- colnames(DM_df)
  initial_overlap <- length(intersect(run_names, rownames(meta_df)))
  if (initial_overlap == 0) {
    pick <- choose_meta_id(meta_df, run_names)
    if (!is.null(pick$best) && pick$score > 0) {
      message("###DEBUG### switching metadata ID column to: ", pick$best,
              " (overlap ", pick$score, ")")
      rownames(meta_df) <- .normalize_runs(meta_df[[pick$best]])
      if (anyDuplicated(rownames(meta_df))) {
        message("###WARNING### duplicates in chosen metadata ID; making them unique")
        rownames(meta_df) <- make.unique(rownames(meta_df))
      }
    }
  }
  
  # ---- Final diagnostics; hard-stop if still zero to avoid empty output ----
  overlap <- length(intersect(colnames(DM_df), rownames(meta_df)))
  message("###DEBUG### data_df dim: ", paste(dim(DM_df), collapse="x"))
  message("###DEBUG### meta_df dim: ", paste(dim(meta_df), collapse="x"))
  message("###DEBUG### overlap (final pre-align): ", overlap)
  
  if (overlap == 0) {
    message("###DEBUG### sample names (data_df, first 5): ",
            paste(utils::head(colnames(DM_df), 5), collapse=", "))
    message("###DEBUG### sample names (meta_df, first 5): ",
            paste(utils::head(rownames(meta_df), 5), collapse=", "))
    stop("No overlapping sample IDs between data columns and metadata rownames after normalization/mapping. ",
         "Check --ids mapping and naming conventions.", call. = FALSE)
  }
  # Optional seed handling for reproducibility
  if (!is.null(params[['seed']])) {
    set.seed(params[['seed']] + as.integer(args$counter))
  }

  # Run with logging
  corrected_df <- run_with_logging(
    function() batch_correct(DM_df, meta_df, params),
    args$log
  )

  # Write output
  corrected_df <- as.data.frame(corrected_df)
  output_file <- file.path(args$out, sprintf("data_%02d.tsv", as.integer(args$counter)))
  fwrite(corrected_df, output_file, row.names = TRUE)
}

if (sys.nframe() == 0) main()
