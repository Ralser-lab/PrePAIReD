###############################################################################
# Normalization module in R
# Provides the following options for normalization: VSN, CyclicLOESS
# Author: Álvaro Gómez Pérez (based on Vadim Farztdinov´s)
# Tested on R version ...
###############################################################################



# Packages
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
  library(vsn)
  library(data.table)
})
source("Scripts/helper_functions_R.R")



#' @description
#' Function to apply variance stabilizing normalization (VSN).   
#' 
#' @param DM_df A \code{data.frame} or numeric matrix, containing log2
#' transformed data. Columns correspond to variables and rows to observations. 
#' @return A \code{data.frame} with the normalized, glog2-transformed data.
#' 
normalize_vsn <- function(DM_df){
  # Basic checks
  if (!is.matrix(DM_df) && !is.data.frame(DM_df)) {
    stop("normalize_vsn(): input must be a matrix or data.frame.", call. = FALSE)
  }
  DM_mat <- as.matrix(DM_df)
  if (!is.numeric(DM_mat)) {
    stop("normalize_vsn(): all values must be numeric.", call. = FALSE)
  }
  
  # Missingness warning (as percentage)
  prop_na <- mean(is.na(DM_mat))
  warning(sprintf(
    "The provided dataset contains %.1f%% missing values. Any value over 20%% will make this normalization method less reliable—consider imputation beforehand or a different normalization method.",
    prop_na * 100
  ))
  
  # Print message informing the user that normalization is starting, with time stamp
  cat(sprintf("\n\n Starting VSN normalization. [%s] \n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # The data provided is log2 transformed, need to be transformed back to 
  # raw intensities for VSN normalization
  DM_mat <- as.matrix(2^DM_mat)
  
  # Perform the normalization
  library(limma)
  DM_mat <- limma::normalizeVSN(DM_mat)
  
  # Print message informing the user that normalization has ended, with time stamp
  cat(sprintf("\n\n Finished VSN normalization. [%s] \n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # Return output
  return(as.data.frame(DM_mat, stringsAsFactors = FALSE))
  
  ### Do we want to keep this verbose possibility?
  # if (normalize$verbose) {
  #   print("Properties of data matrix vsn normalized mean distribution. Summary: ")
  #   print(summary(apply(DM,1,mean, na.rm=TRUE)))
  # }
}



#' @description
#' Apply cyclic LOESS normalization: normalize the columns of a matrix, 
#' cyclically applying LOESS normalization to normalize each pair of columns to 
#' each other.
#' The input DM should contain log2-transformed data.
#' 
#' @param DM_df A \code{data.frame} or numeric matrix, containing log2
#' transformed data. Columns correspond to variables and rows to observations. 
#' @param parameters A list of arguments to pass to \code{limma::normalizeCyclicLoess}:
#' \itemize{
#'    \item \code{span}: The span of LOESS smoothing window, float between 0 and 1. 
#'    Default is \code{NULL}.
#'    \item \code{method_str}: A character string specifying the type of cyclic 
#'    LOESS method to be used. Options are "fast" (default), "affy" and "pairs".
#'    \item \code{iter_int}: An integer, the number of times to cycle through all 
#'    pairs of columns. Default is 3.
#' }
#' 
#' @return A \code{data.frame} with the normalized, log2-transformed data.
normalize_cyclic_loess <- function(DM_df,
                                   parameters) {
  # Basic checks
  if (!is.matrix(DM_df) && !is.data.frame(DM_df)) {
    stop("normalize_vsn(): input must be a matrix or data.frame.", call. = FALSE)
  }
  DM_mat <- as.matrix(DM_df)
  if (!is.numeric(DM_mat)) {
    stop("normalize_vsn(): all values must be numeric.", call. = FALSE)
  }
  
  # Missingness warning (as percentage)
  prop_na <- mean(is.na(DM_mat))
  warning(sprintf(
    "The provided dataset contains %.1f%% missing values. Any value over 20%% will make this normalization method less reliable—consider imputation beforehand or a different normalization method.",
    prop_na * 100
  ))
  
  # Print message informing the user that normalization is starting, with time stamp
  cat(sprintf(
    "\n\n Starting cyclic LOESS normalization (with %s method). [%s] \n",
    parameters[["method_cyclic_loess"]],
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ))
  
  # Make sure the method chosen is a valid option
  method_str <- match.arg(parameters[["method_cyclic_loess"]],
                          choices = c("fast", "affy", "pairs"))
  
  # Perform the normalization
  library(limma)
  DM_df <- limma::normalizeCyclicLoess(DM_mat, 
                                       span = parameters[["span_cyclic_loess"]], 
                                       method = parameters[["method_cyclic_loess"]],
                                       iterations = parameters[["iter_cyclic_loess"]])
  
  # Print message informing the user that normalization has ended, with time stamp
  cat(sprintf("\n\n Finished cyclic LOESS normalization. [%s] \n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # Return output
  return(as.data.frame(DM_mat, stringsAsFactors = FALSE))
}




#' @description
#' Normalize a numeric data matrix or data frame.
#'
#' @param data A numeric data.frame or matrix, log2-transformed.
#'             Columns = samples, rows = precursors.
#' @param method A string: "vsn" or "cyclic_loess".
#' @param params A named list of extra parameters to pass to the method.
#'
#' @return A normalized data.frame with the same dimensions and rownames as `data`.
normalize <- function(data, 
                      method, 
                      params = list()) {
  stopifnot(is.matrix(data) || is.data.frame(data))
  
  # Convert to data.frame for consistent output
  data <- as.data.frame(data)
  
  # Choose method
  main_normalize <- switch(
    tolower(method),
    "vsn"          = function() normalize_vsn(data),
    "cyclic_loess" = function() normalize_cyclic_loess(data, params),
    stop(sprintf("Unknown normalization method: %s", method), call. = FALSE)
  )
  
  # Run normalization (no logging, no file I/O here)
  normalized_df <- main_normalize()
  
  return(normalized_df)
}






#' @title Main function
#' 
#' @description 
#' Parses input from the wrapper, runs normalization and writes output
#' 
#' 
main <- function() {
  # Parse input with functions defined in helper_functions_R.R
  parser <- make_parser()
  args   <- parse_args(parser, args = commandArgs(trailingOnly = TRUE))
  validate_opts(args)
  
  # Read parameters
  parameter_list <- read_JSON(args$params)
  
  # Set seed
  if (!is.null(parameter_list[["seed"]])) {
    set.seed(parameter_list[["seed"]] + args$counter)
  }
  
  # Load data
  DM_df <- as.data.frame(as.matrix(fread(args$data), rownames = 1))
  
  # Call normalize() with logging
  normalized_df <- run_with_logging(
    function() normalize(
      data   = DM_df,
      method = parameter_list[["normalization_method"]],
      params = parameter_list
    ),
    args$log
  )
  
  # Write output
  normalized_df <- as.data.frame(normalized_df)
  output_file <- file.path(args$out, sprintf("data_%d.tsv", as.integer(args$counter)))
  fwrite(normalized_df, output_file, row.names = TRUE)
}

if (sys.nframe() == 0) main()