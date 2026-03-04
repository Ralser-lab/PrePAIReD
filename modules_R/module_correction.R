###############################################################################
# Summarization module in R
# Provides the following options for summarization: medianpolish, plm, fastmaxLFQ
# Author: LCN
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


#' @title Function for peptide-to-protein summarization in R
#' @description
#' Reads a quantitative peptide-level table and a parameter JSON file, then
#' summarizes peptides to protein abundances using one of the supported methods.
#' Wrapped in a logging function to capture console output and messages.
#'
#' @param datapath Path to quantitative input file.
#' @param parampath Path to a JSON parameter file.
#' @param idmappath Path to mapping dataframe.
#' @param outpath Directory where output will be written.
#' @param log_file File path where logs will be written.
#' @param counter Step of the preprocessing pipeline currently run.
summarize <- function(datapath, parampath, outpath, idmappath, log_file, counter){
    
    ## build input paths robustly (no dependence on trailing "/")
    in_data_file <- datapath
    in_mask_file <- sub("\\.tsv$", "_mask.tsv", datapath)
    in_mask_desc_file <- sub("\\.tsv$", "_mask_description.tsv", datapath)
    
    ## read in data
    quant_df <- as.data.frame(read.table(file = in_data_file,
                                         sep = "\t",
                                         header = TRUE,
                                         row.names = 1,
                                         check.names = FALSE))
    
    parameter_list <- read_JSON(parampath)
    
    mapping_df <- as.data.frame(read.table(file = idmappath,
                                           sep = "\t",
                                           header = TRUE,
                                           check.names = FALSE))
    
    ## setting seed
    set.seed(parameter_list[["seed"]] + counter)
    
    ## check if feature and protein columns are in mapping df
    if(!(parameter_list[["feature_col"]] %in% colnames(mapping_df) &
         parameter_list[["protein_col"]] %in% colnames(mapping_df))){
        stop("'feature_col' and/or 'protein_col' provided in parameter list do not exist in id mapping file.")
    }
    
    ## check overlap between row names and mapping
    if(length(intersect(row.names(quant_df),
                        mapping_df[, parameter_list[["feature_col"]]])) == 0){
        stop("No overlap between row names of data and feature ids in id_mapping file.")
    }
    
    ## call summarization function
    main_summarize <- switch(
        tolower(parameter_list[["summarization_method"]]),
        "medianpolish" = function() summarize_median_polish(quant_df, 
                                                            parameter_list,
                                                            mapping_df),
        "fastmaxlfq" = function() summarize_fast_maxlfq(quant_df, 
                                                        parameter_list,
                                                        mapping_df),
        "plm" = function() summarize_plm(quant_df, 
                                         parameter_list,
                                         mapping_df))
    
    protein_df <- run_with_logging(main_summarize, log_file)
    
    ## replace generated NaNs by NAs
    protein_df[is.nan(as.matrix(protein_df))] <- NA
    protein_df[!is.finite(as.matrix(protein_df))] <- NA
    
    ## load mask + mask_description to save unchanged
    mask <- as.data.frame(as.matrix(fread(in_mask_file), rownames = 1))
    mask_description <- as.data.frame(as.matrix(fread(in_mask_desc_file), 
                                                rownames = 1))
    
    ## write outputs
    out_data_file <- file.path(outpath,
                               sprintf("data_%02d.tsv", as.numeric(counter)))
    out_data_file_precursors_unlog <- file.path(outpath,
                                                sprintf("data_%02d_precursors_unlog.tsv",
                                                        as.numeric(counter)))
    out_data_file_unlog <- file.path(outpath,
                                     sprintf("data_%02d_unlog.tsv",
                                             as.numeric(counter)))
    out_mask_file <- file.path(outpath, sprintf("data_%02d_mask.tsv", 
                                                as.numeric(counter)))
    out_mask_desc_file <- file.path(outpath, 
                                    sprintf("data_%02d_mask_description.tsv", 
                                            as.numeric(counter)))
    
    write.table(protein_df,
                out_data_file,
                sep = "\t",
                quote = FALSE,
                col.names = NA)
    
    write.table(2^quant_df,
                out_data_file_precursors_unlog,
                sep = "\t",
                quote = FALSE,
                col.names = NA)
    
    write.table(2^protein_df,
                out_data_file_unlog,
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


#' Summarize precursors/peptides to proteins using median polish
summarize_median_polish <- function(quant_df, parameter_list, mapping_df){
    
    library(iq)
    
    feat_col <- parameter_list[["feature_col"]]
    prot_col <- parameter_list[["protein_col"]]
    
    cat("\n\nStarting median polish summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    quant_list <- split(
        as.data.frame(quant_df),
        mapping_df[match(row.names(quant_df),
                         mapping_df[, feat_col]), prot_col]
    )
    quant_list <- lapply(quant_list, as.matrix)
    
    quant_df_Nfeat <- do.call(rbind, lapply(quant_list, \(x){
        colSums(!is.na(x))
    }))
    
    quant_prot_df <- iq::create_protein_table(quant_list)$estimate
    quant_prot_df[quant_df_Nfeat < parameter_list[["min_feat_per_prot"]]] <- NA
    
    cat("\n\nFinished median polish summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    return(quant_prot_df)
}


#' Summarize precursors/peptides to proteins using fastMaxLFQ
summarize_fast_maxlfq <- function(quant_df, parameter_list, mapping_df){
    
    library(iq)
    
    feat_col <- parameter_list[["feature_col"]]
    prot_col <- parameter_list[["protein_col"]]
    
    cat("\n\nStarting fast MaxLFQ summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    prot_ids <- mapping_df[match(row.names(quant_df), mapping_df[, feat_col]), prot_col]
    
    quant_list <- split(as.data.frame(quant_df), prot_ids)
    
    res <- lapply(names(quant_list), \(p){
        
        x <- as.matrix(quant_list[[p]])
        idx <- which(!is.na(x), arr.ind = TRUE)
        
        if(nrow(idx) == 0) return(NULL)
        
        list(protein = rep(p, nrow(idx)),
             sample = colnames(x)[idx[,2]],
             id = rownames(x)[idx[,1]],
             quant = x[idx])
    })
    
    res <- res[!sapply(res, is.null)]
    
    norm_data <- list(protein_list = unlist(lapply(res, `[[`, "protein"), use.names = FALSE),
                      sample_list = unlist(lapply(res, `[[`, "sample"), use.names = FALSE),
                      id = unlist(lapply(res, `[[`, "id"), use.names = FALSE),
                      quant = unlist(lapply(res, `[[`, "quant"), use.names = FALSE))
    
    quant_prot_df <- iq::fast_MaxLFQ(norm_data)$estimate
    
    cat("\n\nFinished fast MaxLFQ summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    return(quant_prot_df)
}


#' Summarize precursors/peptides to proteins using PLM
summarize_plm <- function(quant_df, parameter_list, mapping_df){
    
    library(preprocessCore)
    
    feat_col <- parameter_list[["feature_col"]]
    prot_col <- parameter_list[["protein_col"]]
    
    cat("\n\nStarting PLM summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    quant_list <- split(
        as.data.frame(quant_df),
        mapping_df[match(row.names(quant_df),
                         mapping_df[, feat_col]), prot_col]
    )
    quant_list <- lapply(quant_list, as.matrix)
    
    quant_df_Nfeat <- do.call(rbind, lapply(quant_list, \(x){
        colSums(!is.na(x))
    }))
    
    sample_names <- colnames(quant_df)
    n_samples <- ncol(quant_df)
    
    quant_prot_df <- do.call(rbind, lapply(quant_list, \(df){
        out <- preprocessCore::rcModelPLM(y = df)
        out$Estimates[seq_len(n_samples)]
    }))
    
    colnames(quant_prot_df) <- sample_names
    quant_prot_df[quant_df_Nfeat < parameter_list[["min_feat_per_prot"]]] <- NA
    
    cat("\n\nFinished PLM summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    return(quant_prot_df)
}


# CLI entrypoint (pipeline integration)
main <- function() {
    parser <- make_parser()
    args <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))
    validate_opts(args)
    
    summarize(datapath = args$data,
              parampath = args$params,
              outpath = args$out,
              idmappath = args$ids,
              log_file = args$log,
              counter = as.numeric(args$counter))
}

if (sys.nframe() == 0) main()