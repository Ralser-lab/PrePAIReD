###############################################################################
# Normalization module in R
# Provides the following options for normalization: VSN, CyclicLOESS
# Author: Álvaro Gómez Pérez (based on Vadim Farztdinov´s)
# Tested on R version ...
###############################################################################



# Packages
library(roxygen2)






###############################################################################
############################## Create functions ###############################
###############################################################################

#' @description
#' Function to apply variance stabilizing normalization (VSN).   
#' 
#' @param DM_df A \code{data.frame} or numeric matrix, containing log2
#' transformed data. Columns correspond to variables and rows to observations. 
#' @return A \code{data.frame} with the normalized, glog2-transformed data.
#' 
normalize_vsn <- function(DM_df){
  # The data provided is log2 transformed, need to be transformed back to 
  # raw intensities for VSN normalization
  DM_df <- 2^DM_df
  
  # Perform the normalization
  library(limma)
  DM_df <- limma::normalizeVSN(DM_df)
  
  # Return output
  return(as.data.frame(DM_df))

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
                                   span = 0.7,
                                   method_str = c("fast", "affy", "pairs"),
                                   iter_int = 3) {
  # Make sure the method chosen is a valid option
  method_str <- match.arg(method_str)
  
  # Perform the normalization
  library(limma)
  DM_df <- limma::normalizeCyclicLoess(DM_df, 
                                span = 0.7, 
                                iterations = iter_int, 
                                method = method_str)
  
  # Return output
  return(as.data.frame(DM_df))
}


#' @description
#' Function for normalizing data in R. Reads the data from the provided location,
#' applies the normalization method specified in the parameter list as "method" 
#' (providing the necessary parameters from this list as well), and writes the
#' resulting, normalized data frame to the specified output location. 
#' 
#' @param data_path A string pointing to the location of the dataset to be normalized
#' This must be a \code{data.frame} or numeric matrix, containing log2 transformed
#' data. Columns correspond to variables and rows to observations. 
#' @param out_dir_path A string pointing to the directory where the produced, 
#' normalized data set should be written as a .tsv file. 
#' @param parameters_path A string pointing to the location of a JSON file, which
#' contains the list of parameters defining how to perform the normalization:
#' \itemize{
#'    \item \code{method} A string defining the normalization method to be used.
#'    Can be "vsn" or "cyclic_loess"
#'    \item Parameters passed to the function normalize_cyclic_loess
#' }
#' 
#' @return Returns \code{NULL}. Writes the normalized version of the dataframe,
#' as a .tsv file, to the location where the original data was located
main_normalize <- function(data_path,
                           out_dir_path,
                           parameters_path) {
  # First of all we make sure the method chosen is one of the options we provide
  method_str <- match.arg(method_str, choices = c("vsn", "cyclic_loess"))
  
  # Read the data 
  DM_df <- as.data.frame(fread(data_path))
  
  # Read the parameters
  library(jsonlite)
  parameter_list <- jsonlite::fromJSON(parameters_path)
  
  # Depending on which method was chosen, perform that normalization
  normalized_df <- switch(
    tolower(parameter_list[["method"]]),
    "vsn" = normalize_vsn(DM_df),
    "cyclic_loess" = normalize_cyclic_loess(DM_df, parameter_list)
  )
  
  # OLD  VERSION - SUBSTITUTED BY THE "SWITCH" ABOVE
  ## Depending on which method was chosen, perform that normalization
  #if (method_str == "vsn") {
  #  DM_df <- normalize_vsn(DM_df)
  #}
  #else if (method_str == "cyclic_loess") {
  #  # If cyclic LOESS was selected, make sure we received the right arguments
  #  cyclic_loess_method_str <- match.arg(parameter_list[["cyclic_loess_method"]])
  #  
  #  # Then perform the normalization
  #  DM_df <- normalize_cyclic_loess(DM_df,
  #                                  span = parameter_list[["span"]],
  #                                  method_str = parameter_list[["method"]],
  #                                  iter_int = parameter_list[["iter"]])
  #}
  
  # Write output to the provided directory
  fwrite(normalized_df, paste(out_dir_path, "dataset_name", "_normalized.tsv", sep=""))
}


















############################ Testing - delete this ############################
# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Prepare data
DM <- as.data.frame(fread("C:/Users/algo12/Downloads/Test_data_reduced.tsv"))
DM_wide <- DM %>%
  dplyr::select(File.Name, Precursor.Normalised, Precursor.Id) %>%
  pivot_wider(values_from = Precursor.Normalised,                               # Probably should do imputation before this, oh well
              names_from = File.Name) %>%
  na.omit() %>%
  as.data.frame()
rownames(DM_wide) <- DM_wide$Precursor.Id
DM_wide <- DM_wide %>%
  dplyr::select(-Precursor.Id) %>%
  log2()

# Normalize VSN
normalized_test_vsn <- main_normalize(DM_df = DM_wide,
                                  method_str = "vsn")

# Normalize cyclic LOESS default parameters
normalized_test_cyclic_loess <- main_normalize(DM_df = DM_wide,
                                               method_str = "cyclic_loess",
                                               cyclic_loess_method_str = "fast",
                                               cyclic_loess_iter_int = 15,
                                               cyclic_loess_weight_vect = rnorm(n = nrow(DM_wide)))


# Compare 
## Get long versions for plotting
DM_wide$Precursor.Id <- rownames(DM_wide)
raw_long <- pivot_longer(DM_wide, 
                         cols = -Precursor.Id,
                         values_to = "Precursor.Normalised",
                         names_to = "Sample")

normalized_test_vsn$Precursor.Id <- rownames(normalized_test_vsn)
norm_vsn_long <- pivot_longer(normalized_test_vsn, 
                              cols = -Precursor.Id,
                              values_to = "Precursor.Normalised",
                              names_to = "Sample")

normalized_test_cyclic_loess$Precursor.Id <- rownames(normalized_test_cyclic_loess)
norm_cyclic_loess_long <- pivot_longer(normalized_test_cyclic_loess, 
                                       cols = -Precursor.Id,
                                       values_to = "Precursor.Normalised",
                                       names_to = "Sample")

## Plot
ggplot(data = raw_long, aes(x = Sample, y = Precursor.Normalised)) +
  geom_boxplot(outliers = F) +
  geom_hline(yintercept = median(raw_long$Precursor.Normalised), col = "orange")
  theme_light() +
  theme(axis.text.x=element_blank()) +
  labs(title = "Raw")
ggplot(data = norm_vsn_long, aes(x = Sample, y = Precursor.Normalised)) +
  geom_boxplot(outliers = F) +
  geom_hline(yintercept = median(norm_vsn_long$Precursor.Normalised), col = "orange") +
  theme_light() +
  theme(axis.text.x=element_blank()) +
  labs(title = "VSN")
ggplot(data = norm_cyclic_loess_long, aes(x = Sample, y = Precursor.Normalised)) +
  geom_boxplot(outliers = F) +
  geom_hline(yintercept = median(norm_cyclic_loess_long$Precursor.Normalised), col = "orange") +
  theme_light() +
  theme(axis.text.x=element_blank()) +
  labs(title = "Cyclic LOESS")















