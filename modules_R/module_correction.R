###############################################################################
# Batch correction module in R
# Provides the following options for batch correction:
# limma, combat_r
# Author: Roza Mizrak @rozsurmiz (pipeline style refactor by LCN)
###############################################################################

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

# Pipeline-safe helper sourcing
script_path <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg)))
    if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(sys.frames()[[1]]$ofile))
    stop("Can't determine script path.")
})()

script_dir <- dirname(script_path)
source(file.path(script_dir, "helper_functions_R.R"))


#' @title Function for batch correction in R
#' @description
#' Reads a quantitative data table and metadata, then performs batch correction
#' using one of the supported methods.
#'
#' The batch correction step is wrapped in a logging function to capture console
#' output and messages to a log file.
#'
#' Additionally, it passes through the corresponding mask and mask_description
#' files unchanged so the pipeline continues to work consistently.
#'
#' @param datapath Path to quantitative input data (.tsv).
#' @param parampath Path to a JSON parameter file.
#' @param outpath Output directory where results are written.
#' @param metadata Path to metadata file.
#' @param log_file Path to log file.
#' @param counter Step of the preprocessing pipeline currently run.
#'
#' @return NULL (writes corrected data and mask files to disk).
batch_correct <- function(datapath, parampath, outpath,
                          metadata, log_file, counter){

    ## build input paths robustly (no dependence on trailing "/")
    in_data_file <- datapath
    in_mask_file <- sub("\\.tsv$", "_mask.tsv", datapath)
    in_mask_desc_file <- sub("\\.tsv$", "_mask_description.tsv", datapath)

    parameter_list <- read_JSON(parampath)

    ## read data (read_tsv_pd checks the feature column for duplicates)
    quant_df <- read_tsv_pd(in_data_file)

    ## read metadata
    ## --- CHANGED: the original code used read.table(row.names = 1), which keys
    ## --- the metadata on its FIRST column (here "Run", the .raw filenames).
    ## --- But the data matrix columns are SAMPLE ids (metadata_sample_column,
    ## --- e.g. "CTT_SID"). Keying on "Run" makes
    ## ---   intersect(colnames(quant_df), rownames(meta_df))
    ## --- empty, so meta_df collapses to 0 rows, the batch factor has 0 levels,
    ## --- and limma dies with "not enough degrees of freedom to define
    ## --- contrasts". We now key the metadata on the configured sample column
    ## --- so it aligns with the data columns (falling back to column 1 if the
    ## --- parameter is absent, preserving the old behaviour for old inputs).
    meta_df <- read_tsv_pd(metadata, row_names = FALSE)

    sample_col <- parameter_list[["metadata_sample_column"]]
    if (is.null(sample_col) || !(sample_col %in% colnames(meta_df))) {
        sample_col <- colnames(meta_df)[1]   # fallback = original behaviour
    }

    check_unique_ids(as.character(meta_df[[sample_col]]), metadata, sample_col)
    rownames(meta_df) <- as.character(meta_df[[sample_col]])

    ## batch handling from schema
    if (!is.null(parameter_list[["batch_cols_meta"]])) {

        batch_cols <- parameter_list[["batch_cols_meta"]]

        if (length(batch_cols) >= 1) {
            parameter_list[["batch"]] <- batch_cols[[1]]
        }

        if (length(batch_cols) >= 2) {
            parameter_list[["batch2"]] <- batch_cols[[2]]
        }
    }

    ## set seed
    set.seed(parameter_list[["seed"]] + counter)

    ## align samples
    common <- intersect(colnames(quant_df), rownames(meta_df))

    ## --- ADDED: fail early with a clear message if data and metadata don't
    ## --- share any samples. Without this, an empty overlap silently produces a
    ## --- 0-column matrix / 0-level batch factor and surfaces later as the
    ## --- cryptic limma "not enough degrees of freedom to define contrasts".
    if (length(common) == 0) {
        stop(sprintf(
            paste0("No overlap between data columns and metadata sample ids ",
                   "(keyed on '%s').\n  Example data column:  '%s'\n  Example ",
                   "metadata id:    '%s'\n  Check that metadata_sample_column ",
                   "matches the sample names used as data columns."),
            sample_col,
            if (ncol(quant_df)) colnames(quant_df)[1] else "<none>",
            if (nrow(meta_df))  rownames(meta_df)[1] else "<none>"),
            call. = FALSE)
    }

    quant_df <- quant_df[, common, drop = FALSE]
    meta_df <- meta_df[common, , drop = FALSE]

    ## method naming fallback
    method <- parameter_list[["correct_method"]]

    ## call batch correction
    main_batch <- switch(
        tolower(method),
        "limma"  = function() batchcorrect_limma(quant_df, meta_df, parameter_list),
        "combat_r" = function() batchcorrect_combat(quant_df, meta_df, parameter_list)
        #"harmony_r"= function() batchcorrect_harmony(quant_df, meta_df, parameter_list) # harmony excluded from the wrapper
    )

    corrected_df <- run_with_logging(main_batch, log_file)

    ## load mask + mask_description to save unchanged
    mask <- read_tsv_pd(in_mask_file)
    mask_description <- read_tsv_pd(in_mask_desc_file, id_label = "filter")

    ## write outputs
    out_data_file <- file.path(outpath,
                               sprintf("data_%02d.tsv", as.numeric(counter)))
    out_mask_file <- file.path(outpath,
                               sprintf("data_%02d_mask.tsv", as.numeric(counter)))
    out_mask_desc_file <- file.path(outpath,
                                    sprintf("data_%02d_mask_description.tsv",
                                            as.numeric(counter)))

    write_tsv_pd(corrected_df, out_data_file)
    write_tsv_pd(mask, out_mask_file)
    write_tsv_pd(mask_description, out_mask_desc_file)

    return(NULL)
}


#' @title Batch correction using limma
#' @description
#' Performs batch correction using \code{limma::removeBatchEffect}.
#'
#' Supports one or two batch factors, optional covariates,
#' and optional design matrix via the \code{Class} parameter.
#'
#' @param quant_df Quantitative data.frame (features × samples).
#' @param meta_df Metadata data.frame (samples × annotations).
#' @param parameter_list Parameter list from JSON file.
#'
#' @return Corrected data.frame.
batchcorrect_limma <- function(quant_df, meta_df, parameter_list){

    library(limma)

    cat("\n\nStarting limma batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")

    batchF1 <- if (!is.null(parameter_list[["batch"]]))
        factor(meta_df[[parameter_list[["batch"]]]]) else NULL

    batchF2 <- if (!is.null(parameter_list[["batch2"]]))
        factor(meta_df[[parameter_list[["batch2"]]]]) else NULL

    covs <- if (!is.null(parameter_list[["covariates"]]))
        as.matrix(meta_df[, parameter_list[["covariates"]], drop = FALSE]) else NULL

    if (is.null(parameter_list[["Class"]])) {

        corrected <- limma::removeBatchEffect(
            as.matrix(quant_df),
            batch = batchF1,
            batch2 = batchF2,
            covariates = covs
        )

    } else {

        Class <- factor(meta_df[[parameter_list[["Class"]]]])
        designM <- model.matrix(~0 + Class)

        corrected <- limma::removeBatchEffect(
            as.matrix(quant_df),
            batch = batchF1,
            batch2 = batchF2,
            covariates = covs,
            design = designM
        )
    }

    cat("\nFinished limma batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")

    as.data.frame(corrected)
}


#' @title Batch correction using ComBat
#' @description
#' Performs batch correction using \code{sva::ComBat}.
#'
#' Optional covariates can be preserved via a design matrix.
#'
#' @param quant_df Quantitative data.frame (features × samples).
#' @param meta_df Metadata data.frame (samples × annotations).
#' @param parameter_list Parameter list from JSON file.
#'
#' @return Corrected data.frame.
batchcorrect_combat <- function(quant_df, meta_df, parameter_list){

    library(sva)

    cat("\n\nStarting ComBat batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")

    mod <- NULL
    if (!is.null(parameter_list[["combat_cols_preserve_meta"]])) {
        mod <- stats::model.matrix(
            ~ .,
            data = meta_df[, parameter_list[["combat_cols_preserve_meta"]], drop = FALSE]
        )
    }

    par_prior <- parameter_list[["par.prior"]]
    if (is.null(par_prior)) par_prior <- TRUE

    corrected <- sva::ComBat(
        dat = as.matrix(quant_df),
        batch = factor(meta_df[[parameter_list[["batch"]]]]),
        mod = mod,
        par.prior = par_prior,
        prior.plots = FALSE
    )

    cat("\nFinished ComBat batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")

    as.data.frame(corrected)
}


#' @title Batch correction using Harmony
#' @description
#' Performs batch correction using PCA followed by
#' \code{harmony::HarmonyMatrix}, and reconstructs the corrected
#' data in expression space while preserving PCA centering.
#'
#' @param quant_df Quantitative data.frame (features × samples).
#' @param meta_df Metadata data.frame (samples × annotations).
#' @param parameter_list Parameter list from JSON file.
#'
#' @return Corrected data.frame.
#batchcorrect_harmony <- function(quant_df, meta_df, parameter_list){
#
#    library(harmony)
#
#    cat("\n\nStarting Harmony batch correction. ",
#        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
#
#    n_pcs <- parameter_list[["n_pcs"]]
#    if (is.null(n_pcs)) n_pcs <- 20
#
#    pca <- prcomp(t(quant_df), center = TRUE, scale. = FALSE)
#    pcs <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]
#
#    vars_use <- unique(na.omit(c(
#        parameter_list[["batch"]],
#        parameter_list[["covariates"]]
#    )))
#
#    harmony_pcs <- harmony::HarmonyMatrix(
#        data_mat = pcs,
#        meta_data = meta_df,
#        vars_use = vars_use,
#        do_pca = FALSE
#    )
#
#    loadings <- pca$rotation[, seq_len(ncol(harmony_pcs)), drop = FALSE]
#
#    centered <- t(quant_df) -
#        matrix(pca$center,
#               nrow = ncol(quant_df),
#               ncol = nrow(quant_df),
#               byrow = TRUE)
#
#    corrected_centered <- harmony_pcs %*% t(loadings)
#
#    corrected <- t(
#        corrected_centered +
#            matrix(pca$center,
#                   nrow = ncol(quant_df),
#                   ncol = nrow(quant_df),
#                   byrow = TRUE)
#    )
#
#    cat("\nFinished Harmony batch correction. ",
#        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
#
#    as.data.frame(corrected)
#}


# CLI entrypoint
main <- function() {
    parser <- make_parser()
    args <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))
    validate_opts(args)

    batch_correct(datapath = args$data,
                  parampath = args$params,
                  outpath = args$out,
                  metadata = args$md,
                  log_file = args$log,
                  counter = as.numeric(args$counter))
}

if (sys.nframe() == 0) main()
