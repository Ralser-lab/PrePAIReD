###############################################################################
# Normalization module in R
# Provides the following options for normalization: VSN, CyclicLOESS
# Author: Álvaro Gómez Pérez (based on Vadim Farztdinov´s)
# Tested on R version ...
###############################################################################



# Packages
library(roxygen2)
library(limma)
library(vsn)






###############################################################################
############################## Create functions ###############################
###############################################################################

#' @description
#' Function to manage input to the normalization function in R
#' Probably from a JSON file, import it and extract the parameters relevant to 
#' normalization? Still need to figure out how this will go exactly
#' 
#' 
handle_json_input <- function(json_file) {
  
}



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
  require(limma)
  DM_df <- normalizeVSN(DM_df);  # What is this semicolon doing, can I get rid of it?
  
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
#'    \item \code{weight_vect}: A numeric vector containing the "probe" weights. 
#'    Default is \code{NULL}.
#'    \item \code{method_str}: A character string specifying the type of cyclic 
#'    LOESS method to be used. Options are "fast" (default), "affy" and "pairs".
#'    \item \code{iter_int}: An integer, the number of times to cycle through all 
#'    pairs of columns. Default is 3.
#' }
#' 
#' @return A \code{data.frame} with the normalized, log2-transformed data.
normalize_cyclic_loess <- function(DM_df,
                                   weight_vect = NULL,
                                   method_str = c("fast", "affy", "pairs"),
                                   iter_int = 3) {
  # Make sure the method chosen is a valid option
  method_str <- match.arg(method_str)
  
  # Perform the normalization
  DM_df <- normalizeCyclicLoess(DM_df, 
                                weights = weight_vect, 
                                span = 0.7, 
                                iterations = iter_int, 
                                method = method_str)
  
  # Return output
  return(as.data.frame(DM_df))
}


#' @description
#' Function for normalizing data in R. Takes the parameters read by 
#' handle_json_input, and calls the appropriate normalization function. 
#' 
#' @param DM_df A \code{data.frame} or numeric matrix, containing log2
#' transformed data. Columns correspond to variables and rows to observations. 
#' @param parameters A list of arguments to pass to \code{limma::normalizeCyclicLoess}:
#' \itemize{
#'    \item \code{weight_vect}: A numeric vector containing the "probe" weights. 
#'    Default is \code{NULL}.
#'    \item \code{method_str}: A character string specifying the type of cyclic 
#'    LOESS method to be used. Options are "fast" (default), "affy" and "pairs".
#'    \item \code{iter_int}: An integer, the number of times to cycle through all 
#'    pairs of columns. Default is 3.
#' }
#' 
#' @return A \code{data.frame} with the normalized, log2-transformed data.
main_normalize <- function(DM_df, 
                           method_str = c("vsn", "cyclic_loess"),
                           cyclic_loess_weight_vect = NULL,
                           cyclic_loess_method_str = c("fast", "affy", "pairs"),
                           cyclic_loess_iter_int = 3) {
  
  # First of all we make sure the method chosen is one of the options we provide
  method_str <- match.arg(method_str)
  
  # Depending on which method was chosen, perform that normalization
  if (method_str == "vsn") {
    DM_df <- normalize_vsn(DM_df)
  }
  else if (method_str == "cyclic_loess") {
    # If cyclic LOESS was selected, make sure we received the right arguments
    cyclic_loess_method_str <- match.arg(cyclic_loess_method_str)
    
    # Then perform the normalization
    DM_df <- normalize_cyclic_loess(DM_df,
                                     cyclic_loess_weight_vect,
                                     cyclic_loess_method_str,
                                     cyclic_loess_iter_int)
  }
  
  # Return output
  return(DM_df)
}








###############################################################################
############################ Perform normalization ############################
###############################################################################













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















