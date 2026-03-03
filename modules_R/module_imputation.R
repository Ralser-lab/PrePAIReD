###############################################################################
# Imputation module in R
# Provides the following options for imputation: missforest, bpca, ppca, ppm, lcg
# Author: LCN (including functions from VF)
# Tested on R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
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
    
    stop("Can't determine script path. Run via: Rscript --vanilla path/to/script.R")
})()

script_dir <- dirname(script_path)

source(file.path(script_dir, "helper_functions_R.R"))


#' @title Function for imputing data in R
#' @description
#' This function reads a quantitative data table and a parameter JSON file,
#' then performs missing value imputation using one of several supported methods.
#' The imputation process is wrapped in a logging function to capture console output
#' and messages to a log file.
#'
#' Additionally, it passes through the corresponding mask and mask_description
#' files unchanged so the pipeline continues to work consistently.
#'
#' @param datapath Path to a directory containing the quantitative data files.
#' @param parampath Path to a JSON parameter file.
#' @param outpath Output directory where results will be written.
#' @param log_file Path to a file where logs will be written.
#' @param counter Step of the preprocessing pipeline currently run.
impute <- function(datapath, parampath, outpath, log_file, counter){
    
    ## build input paths robustly (no dependence on trailing "/")
    in_data_file <- datapath
    in_mask_file <- sub("\\.tsv$", "_mask.tsv", datapath)
    in_mask_desc_file <- sub("\\.tsv$", "_mask_description.tsv", datapath)
    
    ## read in data
    quant_df <- as.data.frame(read.table(file = in_data_file,
                                         sep  = "\t",
                                         header = TRUE,
                                         row.names = 1,
                                         check.names = FALSE))
    
    parameter_list <- read_JSON(parampath)
    
    ## setting seed
    set.seed(parameter_list[["seed"]] + counter)
    
    ## call imputation function
    main_impute <- switch(
        tolower(parameter_list[["imputation_method"]]),
        "missforest" = function() impute_missForest(quant_df, parameter_list, outpath),
        "bpca_R" = function() impute_pca(quant_df, parameter_list, outpath),
        "ppca" = function() impute_pca(quant_df, parameter_list, outpath),
        "ppm" = function() impute_ppm(quant_df, parameter_list, outpath),
        "lcg" = function() impute_MNAR(quant_df, parameter_list)
    )
    
    imputed_df <- run_with_logging(main_impute, log_file)
    
    ## load mask + mask_description
    mask <- as.data.frame(as.matrix(fread(in_mask_file), rownames = 1))
    mask_description <- as.data.frame(as.matrix(fread(in_mask_desc_file), 
                                                rownames = 1))
    
    ## add 0.99 to the mask for every cell that was imputed
    ### positions that were imputed
    imp <- is.na(quant_df)
    
    ### overlapping rows/cols (keeps mask order)
    r <- match(rownames(mask), rownames(imp))
    c <- match(colnames(mask), colnames(imp))
    
    ### valid overlaps
    r_ok <- !is.na(r)
    c_ok <- !is.na(c)
    
    ### add 0.99 for imputed positions
    mask[r_ok, c_ok] <- mask[r_ok, c_ok] +
        0.99 * imp[r[r_ok], c[c_ok], drop = FALSE]
    
    ## write outputs
    out_data_file <- file.path(outpath, sprintf("data_%02d.tsv", 
                                                as.numeric(counter)))
    out_mask_file <- file.path(outpath, sprintf("data_%02d_mask.tsv", 
                                                as.numeric(counter)))
    out_mask_desc_file <- file.path(outpath, 
                                    sprintf("data_%02d_mask_description.tsv", 
                                            as.numeric(counter)))
    
    ## write imputed data.frame into .tsv file
    write.table(imputed_df,
                out_data_file,
                sep = "\t",
                quote = FALSE,
                col.names = NA)
    
    fwrite(mask, out_mask_file, row.names = TRUE, sep = "\t")
    fwrite(mask_description, out_mask_desc_file, row.names = TRUE, sep = "\t")
    
    return(NULL)
}


#' @description
#' Imputation function using the \code{missForest} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{missForest}, such as:
#' \itemize{
#' \item \code{maxiter}: Maximum number of iterations before stopping.
#' \item \code{ntree}: Number of trees to grow in each forest.
#' \item \code{export_output_missForest}: Logical, whether to save the full output object.
#' \item \code{outpath}: Output directory for saving results.
#' }
#'
#' @return A data frame with imputed values.
impute_missForest <- function(quant_df, parameter_list, outpath){
    
    ## load package, throws error if package is not installed
    library(missForest)
    
    ## perform random forest based imputation
    cat("\n\nStarting missforest imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    ## optional: add parallelization, on my WS not faster
    #library(doParallel)
    #cl <- makeCluster(4)
    #registerDoParallel(cl)
    
    imputed_list <- missForest::missForest(
        xmis = quant_df,
        maxiter = parameter_list[["maxiter"]],
        ntree = parameter_list[["ntree"]],
        #parallelize = "forests",
        #maxnodes = 4,
        verbose = TRUE,
        variablewise = parameter_list[["export_output_missForest"]]
    )
    
    #stopCluster(cl)
    
    cat("\nFinished missforest imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    if(parameter_list[["export_output_missForest"]]){
        cat("\n Writing complete output of missForest to output directory.\n")
        saveRDS(imputed_list,
                file.path(outpath,
                       "output_missForest_list.RData"))
    }
    
    imputed_df <- as.data.frame(imputed_list$ximp) # ensures output is a data.frame
    
    return(imputed_df)
}


#' @description
#' Imputation function using the \code{bpca} or \code{ppca} function from the
#' \code{pcaMethods} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{pcaMethods::pca}
#' including:
#' @param parameter_list A list of arguments to pass to \code{pcaMethods::pca}, such as:
#' \itemize{
#' \item \code{imputation_method}: PCA method to use (e.g. \code{"bpca"}, \code{"ppca"}).
#' \item \code{nPCs_imputeR}: Number of principal components for imputation.
#' \item \code{export_output_pca}: Logical indicating whether to save the PCA model object.
#' \item \code{outpath}: Output directory for saving results.
#' }
#'
#' @return A data frame with imputed values.
impute_pca <- function(quant_df, parameter_list, outpath){
    
    ## load package, throws error if package is not installed
    library(pcaMethods)
    
    ## if number of PCs is not defined: sets number of PCs to 15,
    ## if nSamples <15, set nPCs = nSamples
    
    if(is.null(parameter_list[["nPCs_imputeR"]])){
        nPCs <- if (ncol(quant_df) > 15) {
            15
        } else {
            ncol(quant_df)
        }
    } else{
        nPCs <- parameter_list[["nPCs_imputeR"]]
    }
    
    ## perform Bayesian PCA
    cat("\n\nStarting ", parameter_list[["imputation_method"]], " imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    imputed_list <- pcaMethods::pca(
        object = t(quant_df),
        method = parameter_list[["imputation_method"]],
        nPcs = nPCs
    )
    
    cat("\nFinished ", parameter_list[["imputation_method"]], " imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    
    if(parameter_list[["export_output_pca"]]){
        saveRDS(imputed_list, file.path(outpath,
                paste0("output_", parameter_list[["imputation_method"]],"_list.RData")))
    }
    
    ## get the estimated complete observations
    imputed_df <- pcaMethods::completeObs(imputed_list)
    imputed_df <- as.data.frame(t(imputed_df))
    
    return(imputed_df)
}


#' @description
#' Imputation function using \code{ppm} from the \code{mice} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{mice::mice} including:
#' \itemize{
#' \item \code{n_multi_impute}: numeric, number of multiple imputations.
#' }
#'
#' @return A data frame with imputed values.
impute_ppm <- function(quant_df, parameter_list, outpath){
    
    ## load package, throws error if package is not installed
    library(mice)
    #library(parallel)
    
    ## ppm imputation run in the mice package -- slow
    cat("\n\nStarting ppm imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    invisible(capture.output({
        imputed_list <- mice::mice(
            quant_df,
            method = "pmm",
            m = parameter_list[["n_multi_impute"]]
        )
    }))
    
    cat("\nFinished ppm imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    if(parameter_list[["export_output_ppm"]]){
        saveRDS(imputed_list,
                file.path(outpath, "output_ppm_list.RData"))
    }
    
    imputed_df <- as.data.frame(complete(imputed_list))
    
    return(imputed_df)
}


#' @description
#' Imputation function for imputing left-censored missing data by random draws
#' from a Gaussian distribution centered in a minimal value.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass including:
#' \itemize{
#' \item \code{center}: numeric, center of missing values distribution. Default is NULL.
#' \item \code{sigma}: numeric, width on missing values distribution. Default is NULL.
#' \item \code{cutoff_feat_complete}: numeric, detection cutoff. Default is 0.5.
#' }
#'
#' @return A data frame with imputed values.
impute_MNAR <- function(quant_df, parameter_list){
    
    cat("\n\nStarting lcg imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    nSamp <- ncol(quant_df)
    probObserved <- rowSums(!is.na(quant_df))/nSamp
    
    ## get center of missing values distribution
    if (is.null(parameter_list[["center"]])){
        center <- median(apply(quant_df, 2, quantile, prob = 0.01, na.rm = T))
    }
    
    ## get the width on missing values distribution
    if (is.null(parameter_list[["sigma"]])){
        featureSD <- apply(quant_df[probObserved > 0.5, ], 1, sd)
        sigma <- median(featureSD, na.rm = TRUE)
    }
    
    cat("\nCentre of lcg distribution = ", center, "and sigma = ", sigma, ".\n")
    
    ## defining features to be imputed
    impute_subset_df <- quant_df[probObserved < parameter_list[["cutoff_feat_complete"]], ]
    nonimpute_subset_df <- quant_df[probObserved >= parameter_list[["cutoff_feat_complete"]], ]
    
    
    ## imputation of NAs in subset df
    nFeat_subset <- nrow(impute_subset_df)
    impute_subset_df <- apply(impute_subset_df, 2, \(x){
        imputed_vec <- abs(rnorm(nFeat_subset, mean = center, sd = sigma))
        x[is.na(x)] <- imputed_vec[is.na(x)]
        return(x)
    })
    impute_subset_df <- as.data.frame(impute_subset_df)
    
    ## join subsets of data together, sort in original order
    imputed_df <- rbind(impute_subset_df,
                        nonimpute_subset_df)[row.names(quant_df),]
    
    cat("\nFinished lcg imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    return(imputed_df)
}


# CLI entrypoint (pipeline integration)
main <- function() {
    parser <- make_parser()
    args <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))
    validate_opts(args)
    
    impute(datapath = args$data,
           parampath = args$params,
           outpath = args$out,
           log_file = args$log,
           counter = as.numeric(args$counter))
}

if (sys.nframe() == 0) main()
