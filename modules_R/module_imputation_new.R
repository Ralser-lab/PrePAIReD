source("~/Documents/Science/Projects/RalserLabPreprocessing/Preprocessing/modules_R/helper_functions_R.R") ## To-Do: needs to be altered later


#' @title Function for imputing data in R
#' @description
#' This function reads a quantitative data table and a parameter JSON file,
#' then performs missing value imputation using one of several supported methods.
#' The imputation process is wrapped in a logging function to capture console output
#' and messages to a log file.
#'
#' @param datapath Path to a tab-delimited file containing the quantitative 
#' data with missing values. Columns correspond to variables and rows to
#' observations.
#' @param parampath Path to a JSON parameter file. Has to complete a element 
#' /code{method} indicating the imputation method to use. Supported values:
#' `"missforest"`, `"bpca"`, `"ppca"`, `"ppm"`, `"mnar"`.
#' @param log_file Path to a file where logs  will be written.
#' @param counter Step of the preprocessing pupeline currently run
impute <- function(datapath, parampath, outpath, log_file, counter){
    
    ## read in data
    quant_df <- as.data.frame(
        read.table(
            file = paste0(datapath, "data_", as.numeric(counter)-1, ".tsv"),
            sep = '\t',
            header = TRUE)
        )
    parameter_list <- read_JSON(parampath)
    
    ## setting seed
    set.seed(parameter_list[["seed"]] + counter)

    ## call imputation function
    main_impute <- switch(
        tolower(parameter_list[["method_impute"]]),
        "missforest" = function() impute_missForest(quant_df, parameter_list),
        "bpca" = function() impute_pca(quant_df, parameter_list),
        "ppca" = function() impute_pca(quant_df, parameter_list),
        "ppm" = function() impute_ppm(quant_df, parameter_list),
        "mnar" = function() impute_MNAR(quant_df, parameter_list)
    )
    imputed_df <- run_with_logging(main_impute, log_file)
    
    ## write imputed data.frame into .tsv file
    write.table(imputed_df, 
                paste0(outpath, "data_", as.numeric(counter), ".tsv"),
                sep = '\t')
    
    return(NULL)
}


#' @description
#' Imputation function using the \code{missForest} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{missForest}, such as:
#' \itemize{
#'   \item \code{maxiter}: Maximum number of iterations before stopping.
#'   \item \code{ntree}: Number of trees to grow in each forest.
#' }
#'
#' @return A data frame with imputed values.
impute_missForest <- function(quant_df, parameter_list){
    
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
                paste0(parameter_list[["outpath"]], 
                       "output_missForest_list.RData"))
    }
    
    imputed_df <- as.data.frame(imputed_list$ximp)  # ensures output is a data.frame
    
    return(imputed_df)
}


#' @description
#' Imputation function using the \code{bpca} or  \code{ppca} function from the
#' \code{pcaMethods} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{pcaMethods::pca}
#' including:
#' \itemize{
#'   \item \code{}:
#' }
#'
#' @return A data frame with imputed values.
impute_pca <- function(quant_df, parameter_list){
    
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
    cat("\n\nStarting ", parameter_list[["method_impute"]], " imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    imputed_list <- pcaMethods::pca(
        object = t(quant_df),
        method = parameter_list[["method_impute"]],
        nPcs = nPCs
        )
    
    cat("\nFinished ", parameter_list[["method_impute"]], " imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    
    if(parameter_list[["export_output_pca"]]){
        saveRDS(imputed_list,
                paste0(parameter_list[["outpath"]], "output_",
                       parameter_list[["method_impute"]],"_list.RData"))
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
#'   \item \code{n_multi_impute}: numeric, number of multiple imputations.
#' }
#'
#' @return A data frame with imputed values.
impute_ppm <- function(quant_df, parameter_list){
    
    ## load package, throws error if package is not installed
    library(mice)
    #library(parallel)
    
    # ## define numnber of cores
    # ncores_num <- min(parallel::detectCores() - 1, 
    #                   parameter_list[["n_multi_impute"]]) 
    # 
    # if(ncores_num == 0){
    #     ncores_num <- 1
    # }
    # if(ncores_num == 1){
    #     message("Running ppm imputation without parallelization as not enough
    #             cores exist on the machine you are running this script on./n")
    # }
    # else{
    #     message("Running ppm imputation with parallelization.")
    # }
    
    ## ppm imputation run in the mice package -- slow
    cat("\n\nStarting ppm imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n")
    
    ## Option with parall - on my WS did not increase speed, neet to test this
    ## in the final setup and maybe adjust some parameters
    # imputed_list <- mice::futuremice(quant_df,
    #                                  method = "pmm",
    #                                  m = parameter_list[["n_multi_impute"]],
    #                                  parallelseed = parameter_list[["seed"]],
    #                                  n.core = ncores_num)
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
                paste0(parameter_list[["outpath"]], "output_ppm_list.RData"))
    }
    
    imputed_df <- as.data.frame(complete(imputed_list))
    
    return(imputed_df)
}


#' @description
#' Imputation function for imputing left-censored missing data by random draws
#' from a Gaussian distribution centered in a minimal value. Considering a
# peptide/protein expression data matrix with n columns corresponding to
# biological samples and p lines corresponding to peptides/proteins, for each
# sample (column), the mean value of the Gaussian distribution is set to a
# minimal value observed in that sample. The minimal value observed is estimated
# as being the q-th quantile (e.g. q = 0.01) of the observed values in that
# sample. The standard deviation is estimated as the median of the
# peptide/protein-wise standard deviations. Note that when estimating the
# standard deviation of the Gaussian distribution, only the peptides/proteins
# which present more than 50% recorded values are considered.
#
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter_list A list of arguments to pass to \code{mice::mice} including:
#' \itemize{
#'   \item \code{center}: numeric, center of missing values distribution. Default is
#'   \code{NULL}, parameter is then estimated in function.
#'   \item \code{sigma}: numeric, width on missing values distribution. Default is
#'   \code{NULL}, parameter is then estimated in function.
#'   \item \code{cutoff_feat_complete}: numeric, if a feature is detected in less
#'   proportion of samples it will be imputed with MNAR. Features detected in
#'   more samples will not be altered. Default is 0.5. If set to 1, all features
#'   will be imputed
#' }
#'
#' @return A data frame with imputed values.
### imputation function for MNAR
impute_MNAR <- function(quant_df, parameter_list){
    
    cat("\n\nStarting MNAR imputation. ",
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
    
    cat("\nCentre of MNAR distribution = ", center, "and sigma = ", sigma, ".\n")
    
    ## defining features to be imputed
    impute_subset_df <- quant_df[probObserved < parameter_list[["cutoff_feat_complete"]], ]
    nonimpute_subset_df <- quant_df[probObserved >= parameter_list[["cutoff_feat_complete"]], ]
    
    
     #cat("\nFeauteres with less than ", parameter_list[["cutoff_feat_complete"]], 
     #"proportion detected imputed.\nImputing ", nrow(impute_subset_df), 
     #" features from ", nrow(quant_df), " total features.\n")
     
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
    
    cat("\nFinished MNAR imputation. ",
        format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), "\n\n")
    
    return(imputed_df)
}
