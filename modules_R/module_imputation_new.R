source("~/Documents/Science/Projects/RalserLabPreprocessing/helper_functions_R.R")


#' @title Function for imputing data in R
#' @description
#' This function reads a quantitative data table and a parameter JSON file,
#' then performs missing value imputation using one of several supported methods.
#' The imputation process is wrapped in a logging function to capture console output
#' and messages to a log file.
#'
#' @param quant_file Path to a tab-delimited file containing the quantitative 
#' data with missing values. Columns correspond to variables and rows to
#' observations.
#' @param parameter_file Path to a JSON parameter file.
#' @param method A string indicating the imputation method to use.
#'   Supported values: `"missforest"`, `"bpca"`, `"ppca"`, `"ppm"`, `"mnar"`.
#' @param log_file Path to a file where logs  will be written.
impute <- function(quant_file, parameter_file, method, log_file){
    
    ## read in data
    quant_df <- as.data.frame(read.table(file = quant_file,
                                         sep = '\t',
                                         header = TRUE))
    parameter_list <- read_JSON(parameter_file)
    
    ## setting seed
    #set.seed(parameter_list[["seed"]]+parameter_list[["step"]])
    set.seed(42)
    
    ## call imputation function
    main_impute <- switch(tolower(method),
                          "missforest" = function() impute_missForest(quant_df, parameter_list),
                          "bpca" = function() impute_pca(quant_df, parameter_list, method),
                          "ppca" = function() impute_pca(quant_df, parameter_list, method),
                          "ppm" = function() impute_ppm(quant_df, parameter_list),
                          "mnar" = function() impute_MNAR(quant_df, parameter_list))
    imputed_df <- run_with_logging(main_impute, log_file)
    
    ## write imputed data.frame into .tsv file
    #write.table(imputed_df, "filePathAndFileName", sep = '\t')
    
    #return(NULL)
    return(imputed_df)
}


#' @description
#' Imputation function using the \code{missForest} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter A list of arguments to pass to \code{missForest}, such as:
#' \itemize{
#'   \item \code{maxiter}: Maximum number of iterations before stopping.
#'   \item \code{ntree}: Number of trees to grow in each forest.
#' }
#'
#' @return A data frame with imputed values.
impute_missForest <- function(quant_df, parameter){
    
    ## load package, throws error if package is not installed
    library(missForest)
    
    ## perform random forest based imputation
    cat(sprintf("\n\n Starting missForest imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    imputed_list <- missForest::missForest(xmis = quant_df,
                                           maxiter = parameter[["maxiter"]],
                                           ntree = parameter[["ntree"]],
                                           verbose = TRUE) # add parallelization option?
    
    cat(sprintf("\n\n Finished missForest imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    imputed_df <- as.data.frame(imputed_list$ximp)  # ensures output is a data.frame
    
    return(imputed_df)
}
#' @description
#' Imputation function using the \code{bpca} or  \code{ppca} function from the
#' \code{pcaMethods} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter A list of arguments to pass to \code{pcaMethods::pca}
#' including:
#' \itemize{
#'   \item \code{}:
#' }
#' @param method string, defining pcaMethod to be used. Has to be defined as
#' 'bpca' or 'ppca'
#'
#' @return A data frame with imputed values.
impute_pca <- function(quant_df, parameter, method){
    
    ## load package, throws error if package is not installed
    library(pcaMethods)
    
    ## if number of PCs is not defined: sets number of PCs to 15,
    ## if nSamples <15, set nPCs = nSamples
    if(is.null(parameter[["nPCs_imputeR"]])){
        nPCs <- if (ncol(mat) > 15) {
            15
        } else {
            ncol(mat)
        }
    } else{
        nPCs <- parameter[["nPCs_imputeR"]]
    }
    
    ## perform Bayesian PCA
    cat(sprintf("\n\n Starting ", method, " imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    imputed_list <- pcaMethods::pca(xmis = t(quant_df),
                                    method = method,
                                    nPcs = nPCs)
    
    cat(sprintf("\n\n Finished ", method, " imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    ## get the estimated complete observations
    imputed_df <- mice::completeObs(imputed_list)
    imputed_df <- as.data.frame(t(imputed_df))
    
    return(imputed_df)
}
#' @description
#' Imputation function using \code{ppm} from the \code{mice} package.
#'
#' @param quant_df A \code{data.frame} or numeric matrix with missing values.
#' Columns correspond to variables and rows to observations.
#' @param parameter A list of arguments to pass to \code{mice::mice} including:
#' \itemize{
#'   \item \code{n_multi_impute}: numeric, number of multiple imputations.
#' }
#'
#' @return A data frame with imputed values.
impute_ppm <- function(quant_df, parameter){
    
    ## load package, throws error if package is not installed
    library(mice)
    
    ## ppm imputation run in the mice package -- slow
    cat(sprintf("\n\n Starting ppm imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    imputed_list <- mice::mice(quant_df,
                               method = "pmm",
                               m = parameter[["n_multi_impute"]])
    
    cat(sprintf("\n\n Finished ppm imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
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
#' @param parameter A list of arguments to pass to \code{mice::mice} including:
#' \itemize{
#'   \item \code{center}: numeric, center of missing values distribution. Default is
#'   \code{NULL}, parameter is then estimated in function.
#'   \item \code{sigma}: numeric, width on missing values distribution. Default is
#'   \code{NULL}, parameter is then estimated in function.
#'   \item \code{proportionNA}: numeric, if a feature is detected in less
#'   proportion of samples it will be imputed with MNAR. Features detected in
#'   more samples will not be altered. Default is 0.5. If set to 1, all features
#'   will be imputed
#' }
#'
#' @return A data frame with imputed values.
### imputation function for MNAR
impute_MNAR <- function(quant_df, parameter){
    
    cat(sprintf("\n\n Starting missForest imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    nSamp <- ncol(quant_df)
    probObserved <- rowSums(!is.na(quant_df))/nSamp
    
    ## get center of missing values distribution
    if (is.null(parameter[["center"]])){
        center <- median(apply(quant_df, 2, quantile, prob = 0.01, na.rm = T))
    }
    
    ## get the width on missing values distribution
    if (is.null(parameter[["sigma"]])){
        featureSD <- apply(quant_df[probObserved > 0.5, ], 1, sd)
        sigma <- median(featureSD, na.rm = TRUE)
    }
    
    ## defining features to be imputed
    impute_subset_df <- quant_df[probObserved < parameter[["proportionNA"]], ]
    nonimpute_subset_df <- quant_df[probObserved >= parameter[["proportionNA"]], ]
    
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
    
    cat(sprintf("\n\n Finished MNAR imputation. [%s] \n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    
    return(imputed_df)
}