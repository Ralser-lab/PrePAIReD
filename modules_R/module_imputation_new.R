library(jsonlite)
library(mice)
library(missForest)
library(pcaMethods)


#' @description
#' Function for imputing data in R
impute <- function(quant_file, parameter_file, method){
    
    ## read in data
    quant_df <- as.data.frame(read.table(file = quantfile, 
                                         sep = '\t', 
                                         header = TRUE))
    parameter_list <- fromJSON(parameter)
    
    ## call imputation function
    imputed_df <- switch(
        tolower(method),
        "missforest" = impute_missForest(quant_df, parameter_list),
        "bpca" = impute_pca(quant_df, parameter_list, method),
        "ppca" = impute_pca(quant_df, parameter_list, method),
        "ppm" = impute_ppm(quant_df, parameter_list),
        "MNAR" = impute_MNAR(quant_df, parameter_list)
        )
    
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
    
    imputed_list <- missForest::missForest(
        xmis = quant_df, 
        maxiter = parameter[["maxiter"]],
        ntree = parameter[["ntree"]] # add parallelization option?
    )
    
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
    
    ## set number of PCs to 15, if nSamples <15, set nPCs = nSamples
    nPCs <- if (ncol(mat) > 15) {
        15
    } else {
        ncol(mat)
    } 

    ## perform Bayesian PCA 
    imputed_list <- pcaMethods::pca(
        xmis = t(quant_df), 
        method = method,
        nPcs = nPCs
        )
    
    ## get the estimated complete observations
    imputed_df <- completeObs(imputed_list)
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

    ## ppm imputation run in the mice package
    imputed_list <- mice::mice(quant_df, 
                               method = "pmm",
                               m = parameter[["n_multi_impute"]])
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
    
    return(imputed_df)
}
