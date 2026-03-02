###############################################################################
# Batch correction module in R
# Provides the following options for batch correction:
# limma_r, combat_r, harmony_r
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
    
    ## read data
    quant_df <- as.data.frame(read.table(
        file = in_data_file,
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    ))
    
    meta_df <- as.data.frame(read.table(
        file = metadata,
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    ))
    
    parameter_list <- read_JSON(parampath)
    
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
    quant_df <- quant_df[, common, drop = FALSE]
    meta_df <- meta_df[common, , drop = FALSE]
    
    ## method naming fallback
    method <- parameter_list[["correct_method"]]
    
    ## call batch correction
    main_batch <- switch(
        tolower(method),
        "limma_r"  = function() batchcorrect_limma(quant_df, meta_df, parameter_list),
        "combat_r" = function() batchcorrect_combat(quant_df, meta_df, parameter_list),
        "harmony_r"= function() batchcorrect_harmony(quant_df, meta_df, parameter_list)
    )
    
    corrected_df <- run_with_logging(main_batch, log_file)
    
    ## load mask + mask_description to save unchanged
    mask <- as.data.frame(as.matrix(fread(in_mask_file), rownames = 1))
    mask_description <- as.data.frame(as.matrix(fread(in_mask_desc_file),
                                                rownames = 1))
    
    ## write outputs
    out_data_file <- file.path(outpath,
                               sprintf("data_%02d.tsv", as.numeric(counter)))
    out_mask_file <- file.path(outpath,
                               sprintf("data_%02d_mask.tsv", as.numeric(counter)))
    out_mask_desc_file <- file.path(outpath,
                                    sprintf("data_%02d_mask_description.tsv",
                                            as.numeric(counter)))
    
    write.table(corrected_df,
                out_data_file,
                sep = "\t",
                quote = FALSE,
                col.names = NA)
    
    fwrite(mask,
           out_mask_file,
           row.names = TRUE,
           sep = "\t")
    
    fwrite(mask_description,
           out_mask_desc_file,
           row.names = TRUE,
           sep = "\t")
    
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
    if (!is.null(parameter_list[["covariates_preserve"]])) {
        mod <- stats::model.matrix(
            ~ .,
            data = meta_df[, parameter_list[["covariates_preserve"]], drop = FALSE]
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
batchcorrect_harmony <- function(quant_df, meta_df, parameter_list){
    
    library(harmony)
    
    cat("\n\nStarting Harmony batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    n_pcs <- parameter_list[["n_pcs"]]
    if (is.null(n_pcs)) n_pcs <- 20
    
    pca <- prcomp(t(quant_df), center = TRUE, scale. = FALSE)
    pcs <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]
    
    vars_use <- unique(na.omit(c(
        parameter_list[["batch"]],
        parameter_list[["covariates"]]
    )))
    
    harmony_pcs <- harmony::HarmonyMatrix(
        data_mat = pcs,
        meta_data = meta_df,
        vars_use = vars_use,
        do_pca = FALSE
    )
    
    loadings <- pca$rotation[, seq_len(ncol(harmony_pcs)), drop = FALSE]
    
    centered <- t(quant_df) -
        matrix(pca$center,
               nrow = ncol(quant_df),
               ncol = nrow(quant_df),
               byrow = TRUE)
    
    corrected_centered <- harmony_pcs %*% t(loadings)
    
    corrected <- t(
        corrected_centered +
            matrix(pca$center,
                   nrow = ncol(quant_df),
                   ncol = nrow(quant_df),
                   byrow = TRUE)
    )
    
    cat("\nFinished Harmony batch correction. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    as.data.frame(corrected)
}


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