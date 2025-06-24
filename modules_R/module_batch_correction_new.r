#===========================================================
# Unified Batch‐Correction Pipeline in R
#===========================================================

# Install/load required packages
if (!requireNamespace("limma", quietly=TRUE))   install.packages("limma")
if (!requireNamespace("sva", quietly=TRUE))     install.packages("sva")
if (!requireNamespace("harmony", quietly=TRUE)) install.packages("harmony")

library(limma)
library(sva)
library(harmony)

#------------------------------------------------------------------------------  
# 1) batch_effect_removal: LIMMA & ComBat
#------------------------------------------------------------------------------  
#' Batch‐effect removal for LIMMA and ComBat
#'
#' @param Data   list with slots:
#'               $Matrix   numeric matrix (features × samples)
#'               $MetaData data.frame (rows = samples, cols = metadata vars)
#' @param batchER list with settings:
#'               method     "LIMMA" or "ComBat"
#'               batch      metadata column name (string)
#'               batch2     optional second batch factor (string)
#'               covariates optional character vector of metadata cols
#'               Class      optional metadata column to preserve
#'
#' @return Data   same list with $Matrix replaced by corrected matrix
batch_effect_removal <- function(Data, batchER) {
  DM <- Data$Matrix
  md <- Data$MetaData
  method <- toupper(batchER[["method"]])

  if (method == "LIMMA") {
    # extract factors/covariates
    batchF1 <- if (!is.null(batchER$batch)) factor(md[[ batchER$batch ]]) else NULL
    batchF2 <- if (!is.null(batchER$batch2)) factor(md[[ batchER$batch2 ]]) else NULL
    covs    <- if (!is.null(batchER$covariates)) 
                 as.matrix(md[, batchER$covariates, drop=FALSE]) 
               else NULL

    if (is.null(batchER$Class)) {
      DM <- limma::removeBatchEffect(
              DM,
              batch      = batchF1,
              batch2     = batchF2,
              covariates = covs
            )
    } else {
      Class   <- factor(md[[ batchER$Class ]])
      designM <- model.matrix(~0 + Class)
      DM <- limma::removeBatchEffect(
              DM,
              batch      = batchF1,
              batch2     = batchF2,
              covariates = covs,
              design     = designM
            )
    }

  } else if (method == "COMBAT") {
    DM <- sva::ComBat(
            dat       = DM,
            batch     = md[[ batchER$batch ]],
            mod       = NULL,
            par.prior = TRUE,
            prior.plots = FALSE
          )

  } else {
    stop("Unsupported method in batch_effect_removal: ", batchER[["method"]])
  }

  # preserve dims & names
  Data$Matrix <- DM
  return(Data)
}

#------------------------------------------------------------------------------  
# 2) batch_effect_harmony: Harmony integration
#------------------------------------------------------------------------------  
#' Batch‐effect removal via Harmony (PCA‐level integration + back‐projection)
#'
#' @param Data   same as above
#' @param batchER same list; must include:
#'               method = "Harmony"
#'               batch  metadata column name
#'               covariates optional character vector
#'               n_pcs  integer (default 20)
#'               any additional Harmony params (e.g. theta)
#'
#' @return Data   with $Matrix replaced by batch‐corrected matrix
batch_effect_harmony <- function(Data, batchER) {
  DM <- Data$Matrix
  md <- Data$MetaData

  batch_col  <- batchER[["batch"]]
  covar_cols <- batchER[["covariates"]]
  n_pcs      <- batchER[["n_pcs"]] %||% 20

  # 1) PCA
  pca <- prcomp(t(DM), center = TRUE, scale. = FALSE)
  pcs <- pca$x[, seq_len(n_pcs), drop=FALSE]

  # 2) Harmony on PCs
  vars <- c(batch_col, covar_cols)
  # collect extra args for HarmonyMatrix
  extras <- batchER[ setdiff(names(batchER), 
              c("method","batch","batch2","covariates","Class","n_pcs")) ]
  harmony_pcs <- do.call(
    harmony::HarmonyMatrix,
    c(
      list(
        data_mat  = pcs,
        meta_data = md,
        vars_use  = vars,
        do_pca    = FALSE
      ),
      extras
    )
  )

  # 3) back‐project
  loadings <- pca$rotation[, seq_len(n_pcs), drop=FALSE]
  centered <- t(DM) - matrix(pca$center,
                             nrow = ncol(DM),
                             ncol = nrow(DM),
                             byrow = TRUE)
  corrected_centered <- harmony_pcs %*% t(loadings)
  corrected <- t(corrected_centered + 
                 matrix(pca$center,
                        nrow = ncol(DM),
                        ncol = nrow(DM),
                        byrow = TRUE))

  Data$Matrix <- corrected
  return(Data)
}

#------------------------------------------------------------------------------  
# 3) batch_correct: main dispatcher
#------------------------------------------------------------------------------  
#' Unified batch‐correction interface
#'
#' @param Data    list with $Matrix and $MetaData
#' @param batchER list with settings (see above)
#'
#' @return Data    with batch‐corrected $Matrix
batch_correct <- function(Data, batchER) {
  method <- toupper(batchER[["method"]])

  if (method %in% c("LIMMA", "COMBAT")) {
    batch_effect_removal(Data, batchER)
  } else if (method == "HARMONY") {
    batch_effect_harmony(Data, batchER)
  } else {
    stop("Unknown batch correction method: ", batchER[["method"]])
  }
}

#===========================================================
# Usage example
#===========================================================
# Data <- list(
#   Matrix   = my_matrix,        # features × samples
#   MetaData = my_sample_meta,   # samples × metadata vars
#   PepDesc  = my_peptide_desc   # optional
# )
#
# batchER <- list(
#   method     = "LIMMA",       # or "ComBat", "Harmony"
#   batch      =  "Batch",
#   batch2     = NULL,
#   covariates = c("injection_order","operator"),
#   Class      = NULL,
#   n_pcs      = 30,
#   theta      = 2                # harmony hyperparam
# )
#
# Data_corrected <- batch_correct(Data, batchER)
