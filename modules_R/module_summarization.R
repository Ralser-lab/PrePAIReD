source("~/Documents/Science/Projects/RalserLabPreprocessing/Preprocessing/modules_R/helper_functions_R.R") ## To-Do: needs to be altered later

#' @title Function for peptide-to-protein summarization in R
#' @description
#' Reads a quantitative peptide-level table and a parameter JSON file, then
#' summarizes peptides to protein abundances using one of the supported methods.
#' Wrapped in a logging function to capture console output and messages.
#'
#' @param datapath Path to a tab-delimited file containing the quantitative 
#' data with precursor/peptide values that should be summarized to protein 
#' abundances. Columns correspond to variables and rows to observations.
#' @param parampath Path to a JSON parameter file. Has to complete a element 
#' /code{method_summarize} indicating the imputation method to use. Supported values:
#' `"medianpolish"`, `"plm"`.
#' @param idmappath Path to mapping dataframe
#' @param outpath Directory where "data_{counter}.tsv" will be written.
#' @param log_file File path where logs will be written.
#' @param counter Step of the preprocessing pipeline currently run.source("~/Documents/Science/Projects/RalserLabPreprocessing/Preprocessing/modules_R/helper_functions_R.R") ## To-Do: needs to be altered later

summarize <- function(datapath, parampath, outpath, idmappath, og_file, counter){
    
    ## read in data
    quant_df <- as.data.frame(
        read.table(
            file = paste0(datapath, "data_", as.numeric(counter)-1, ".tsv"),
            sep = '\t',
            header = TRUE)
    )
    parameter_list <- read_JSON
    
    mapping_df <- as.data.frame(
        read.table(
            file = idmappath,
            sep = '\t',
            header = TRUE)
    )
    
    ## setting seed
    set.seed(parameter_list[["seed"]] + counter)
    
    ## check if feature and protein columns are in mapping df
    if(!(parameter_list[["feature_col"]] %in% colnames(mapping_df)&
         parameter_list[["protein_col"]] %in% colnames(mapping_df))){
        stop("'feature_col' and/or 'protein_col' provided in parameter list do not
         exist in id mapping file.")
    }
    
    ## check overlap between row names of quant df and feature ids in parameter list
    if(length(intersect(row.names(quant_df), 
                        mapping_df[parameter_list[["feature_col"]],]))==0){
        stop("No overlap between the row names of the data frame in 'datapath'
             and feature ids in the 'feature_col' column of the id_mapping data
             frame.")
    }
    
    ## call summarization function
    main_summarize <- switch(
        tolower(parameter_list[["method_summarize"]]),
        "summarize_median_polish" = function() summarize_medianPolish(quant_df, parameter_list,
                                                           mapping_df),
        "summarize_plm"          = function() summarize_PLM(quant_df, parameter_list,
                                                  mapping_df)
    )
    protein_df <- run_with_logging(main_summarize, log_file)
    
    ## replace generated NaNs by NAs / non-finite by NA (as in reference code)
    protein_df[is.nan(as.matrix(protein_df))] <- NA
    protein_df[!is.finite(as.matrix(protein_df))] <- NA
    
    ## write summarized data.frame into .tsv file
    write.table(protein_df, 
                paste0(outpath, "data_", as.numeric(counter), ".tsv"),
                sep = '\t')
    
    return(NULL)
}


#' Summarize precursors/peptides to proteins using median polish
#'
#' Summarization data to protein abundances using \code{stats::medpolish}.
#'
#' @param quant_df data.frame with peptide-level quantities. Must contain the
#'   feature ID column named \code{parameter_list[["feature"]]} and one or more
#'   numeric sample columns.
#' @param parameter_list list with at least:
#'   \itemize{
#'     \item \code{feature_col} — name of the peptide/precursor column in both
#'       \code{quant_df} and \code{mapping_df}.
#'     \item \code{proteins_col} — name of the protein ID column in \code{mapping_df}.
#'   }
#' @param mapping_df data.frame mapping features (peptides/precursors) to proteins
#' 
#' @return A data.frame of protein (rows) by sample (columns) abundances.
summarize_median_polish <- function(quant_df, parameter_list, mapping_df){
    
    library(iq)
    
    feat_col <- parameter_list[["feature_col"]]
    prot_col <- parameter_list[["protein_col"]]
    
    cat("\n\nStarting median polish summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    ## split data.frame into one matrice per protein
    quant_list <- split(as.data.frame(quant_df), 
                        mapping_df[match(row.names(quant_df), 
                                         mapping_df[, feat_col]), prot_col])
    quant_list <- lapply(quant_list, as.matrix)
    
    ## matrix with number of precursors used for sample
    quant_df_Nfeat <- do.call(rbind, lapply(quant_list, \(x){
        colSums(!is.na(x))
    }))
    
    ## estimate protein abundances & filter for min precursor
    quant_prot_df <- iq::create_protein_table(quant_list, )$estimate
    quant_prot_df[quant_df_Nfeat<parameter_list[["min_feat_per_prot"]]] <- NA
    
    cat("\n\nFinished median polish summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    return(quant_prot_df)
}



#' Summarize precursors/peptides to proteins using plm
#'
#' Summarization data to protein abundances using \code{stats::medpolish}.
#'
#' @param quant_df data.frame with peptide-level quantities. Must contain the
#'   feature ID column named \code{parameter_list[["feature"]]} and one or more
#'   numeric sample columns.
#' @param parameter_list list with at least:
#'   \itemize{
#'     \item \code{feature_col} — name of the peptide/precursor column in both
#'       \code{quant_df} and \code{mapping_df}.
#'     \item \code{protein_col} — name of the protein ID column in \code{mapping_df}.
#'   }
#' @param mapping_df data.frame mapping features (peptides/precursors) to proteins
#' 
#' @return A data.frame of protein (rows) by sample (columns) abundances.
summarize_plm <- function(quant_df, parameter_list, mapping_df){
    
    library(preprocessCore)
    
    feat_col <- parameter_list[["feature_col"]]
    prot_col <- parameter_list[["protein_col"]]
    
    cat("\n\nStarting PLM summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    ## split data.frame into one matrice per protein
    quant_list <- split(as.data.frame(quant_df), 
                        mapping_df[match(row.names(quant_df), 
                                         mapping_df[, feat_col]), prot_col])
    quant_list <- lapply(quant_list, as.matrix)
    
    ## matrix with number of precursors used for sample
    quant_df_Nfeat <- do.call(rbind, lapply(quant_list, \(x){
        colSums(!is.na(x))
    }))
    
    ## estimate protein abundances & filter for min precursor
    sample_names <- colnames(quant_df)
    n_samples <- ncol(quant_df)
    quant_prot_df <- do.call(rbind, lapply(quant_list, \(df){
        out <- preprocessCore::rcModelPLM(y = df)
        out$Estimates[seq(1, n_samples, 1)]
    }))
    colnames(quant_prot_df) <- sample_names
    quant_prot_df[quant_df_Nfeat<parameter_list[["min_feat_per_prot"]]] <- NA
    
    cat("\n\nFinished PLM summarization. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    return(quant_prot_df)
}


