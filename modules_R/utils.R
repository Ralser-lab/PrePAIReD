## Utility functions used to support functions to analyze and visualise ms-dap objects
##
## Program Author: Orr Shomroni, Charité, Core Facility - HTMS, Berlin
##
## Date Created: 2025-02-05
##
## Email: orr.shomroni@charite.de

#' Determine which proteins are expressed over a certain threshold
#'
#' converts a table with peptide-level abundance values to protein-level data by combining data from respective peptides per protein (per sample).
#'
#' @param df Data.frame of peptides expressed in different samples
#' @param exprs_thresh Threshold of expression
exprs_prots <- function(df, exprs_thresh = 0.66) {
  mat <- dfToMat(df, protein_id ~ sample_id, "intensity")
  # Count how many peptides are not NA
  vec <- apply(!is.na(mat), 1, sum) / ncol(mat) > exprs_thresh
  return(vec)
}

#' Circumvent cache_filtering_data function
#' @description
#' The function "cache_filtering_data" from the msdap package creates objects "dt_pep_group" and "groups" in the dataset list, but when the number of samples is large (e.g. ~2000), the "pre-cache CoV" part of the function takes a really long time (CoV is calculated by creating a peptide intensity matrix for all peptides in samples, the matrix is normalized with "vwmb" algorithm and used to calculate the CoV; the "vwmb" normalization takes too long). The function here traverses this normalization and CoV calculation (can be calculated later) so when normalizing the dataset, these objects in the list are there and do not need to be calculated
#' @param dataset msdap list object
#' @returns msdap list object with "dt_pep_group" and "groups"
cache_filtering_data_noCoV <- function(dataset) {
  # add keys to samples table
  s <- data.table::setDT(dataset$samples) # convert to data.table, editing this object BY REFERENCE downstream
  s[, key_sample := .GRP, by = sample_id]
  s[, key_group := .GRP, by = group]
  rm(s)

  # add keys to peptides table
  i <- data.table::chmatch(dataset$peptides$sample_id, dataset$samples$sample_id)
  dataset$peptides$key_sample <- dataset$samples$key_sample[i]
  dataset$peptides$key_group <- dataset$samples$key_group[i]
  rm(i)

  # new reference table for groups
  dataset$groups <- dataset$samples %>%
    group_by(group) %>%
    summarise(size = n(), size_noexclude = sum(!exclude), key_group = key_group[1]) %>%
    arrange(key_group)

  # use a subset of peptides where need to make unique keys
  p <- data.table::setDT(dataset$peptides) # convert to data.table, editing this object BY REFERENCE downstream
  p[, key_peptide := .GRP, by = peptide_id]
  p[, key_protein := .GRP, by = protein_id]
  p[, key_peptide_sample := .GRP, by = list(key_peptide, key_sample)]
  p[, key_protein_group := .GRP, by = list(key_protein, key_group)]
  p[, key_peptide_group := .GRP, by = list(key_peptide, key_group)]
  # flag 'exclude' samples
  kpg <- p$key_peptide_group
  kpg[p$key_sample %in% (dataset$samples %>% filter(exclude == T) %>% pull(key_sample))] <- NA
  p[, key_peptide_group_noexclude := kpg]
  # detect not as a boolean, but scaled by number of detect per sample
  # ` * !is.na(key_peptide_group_noexclude)` is simply to remove all data points for samples that are excluded, which allows us to downstream group by key_group and not worry about excluded samples
  p[, `:=`(
    detect_scaled_noexclude = (detect / sum(detect)) * !is.na(key_peptide_group_noexclude),
    quant_scaled_noexclude = 1 / .N * !is.na(key_peptide_group_noexclude)
  ),
  by = key_sample
  ]

  ## efficiently collapse filtering data for each peptide * 'sample group'
  dataset$dt_pep_group <- p[, .(
    key_peptide = key_peptide[1], key_protein = key_protein[1], key_group = key_group[1], key_protein_group = key_protein_group[1],
    nquant = .N,
    ndetect = sum(detect),
    nquant_noexclude = sum(!is.na(key_peptide_group_noexclude)),
    ndetect_noexclude = sum(!is.na(key_peptide_group_noexclude) & detect),
    ndetect_scaled_noexclude = sum(detect_scaled_noexclude),
    nquant_scaled_noexclude = sum(quant_scaled_noexclude)
  ),
  by = key_peptide_group
  ]
  rm(p)

  # add group sizes and cache the 'fraction of samples within each group' where a peptide is detected
  data.table::setkey(dataset$dt_pep_group, key_group) # set key for data.table::merge()
  dt_criteria <- data.table::data.table(dataset$groups %>% dplyr::select(key_group, size, size_noexclude) %>% mutate_all(as.integer), key = "key_group")
  dataset$dt_pep_group <- data.table::merge.data.table(dataset$dt_pep_group, dt_criteria, all.x = T, sort = FALSE)
  dataset$dt_pep_group[, `:=`(ndetect_fraction = ndetect / size, ndetect_noexclude_fraction = ndetect_noexclude / size_noexclude)]

  dataset$peptides <- as_tibble(dataset$peptides)
  dataset$samples <- as_tibble(dataset$samples)
  dataset$proteins <- as_tibble(dataset$proteins)
  return(dataset)
}

#' @description
#' Remove target peptides from dataset
#' @param dataset msdap dataset
#' @param exclude peptide IDs to be excluded from dataset
#' @return Datasets with target peptide IDs filtered out
remove_exclude_peps <- function(dataset, exclude) {
  dataset$peptides <- dataset$peptides[!(dataset$peptides$peptide_id %in% exclude), ]
  dataset$proteins <- dataset$proteins[dataset$proteins$protein_id %in% unique(dataset$peptides$protein_id), ]
  return(dataset)
}

#' Remove target samples from dataset
#' @param dataset msdap dataset
#' @param exclude sample IDs to be excluded from the dataset
remove_exclude_samps <- function(dataset, exclude) {
  dataset$samples <- dataset$samples[!(dataset$samples$sample_id %in% exclude), ]
  dataset$peptides <- dataset$peptides[dataset$peptides$sample_id %in% dataset$samples$sample_id, ]
  return(dataset)
}

#' @description
#' Impute missing values in a matrix
#'
#' @param mat Matrix with missing values (NAs)
#' @param method Which method to use for imputation
#' @returns imputed matrix
impute.matrix <- function(mat, method = c("rf", "knn", "bpca", "ppca", "pmm")) {
  if (method == "rf") {
    require(missForest)
    mat.impute <- missForest(as.matrix(mat))
    imputed <- mat.impute$ximp
  }
  if (method == "knn") {
    ## Impute method
    require(impute)
    x <- as.matrix(mat)
    colmax <- 0.8
    p <- nrow(x)
    col.nas <- drop(rep(1, p) %*% is.na(x))
    x <- x[, col.nas <= colmax * p]
    mat.impute <- impute.knn(x, colmax = colmax)
    imputed <- mat.impute$data
    ## VIM method
    # imputed <- VIM::kNN(as.matrix(mat), imp_var = FALSE)
  }
  if (method == "bpca" | method == "ppca") {
    require(pcaMethods)
    ## Perform Bayesian PCA with 15 components (unless fewer samples, in which case use the number of samples)
    pc <- pca(t(mat), method = method, nPcs = if (ncol(mat) > 15) {
      15
    } else {
      ncol(mat)
    }) # , seed = 123, maxIterations = 2000)
    
    ## Get the estimated complete observations
    cObs <- completeObs(pc)
    imputed <- t(cObs)
  }
  if (method == "ppm") {
    imputed <- complete(mice::mice(mat, meth = "pmm"))
  }
  return(imputed)
}

#' Remove batch effect
#' @param dataset
#' @param col Column used to extract intensity values from
#' @param batch Factor used to correct the batch effect for (batch, plate or clinical group)
remove_batch_effect <- function(dataset, col, batch = "ms_batch") {
  require(limma)
  # Convert data.frame of intensities for peptides in samples to expression matrix
  mat <- dfToMat(dataset$peptides, col = col)
  
  # Remove batch effect
  # NOTE: for the time being, use limma to remove batch effect, but in the long run it should not be used because it uses a linear model (lmFit function) to remove the batch effect, but it is not ideal for proteomic data which often has missing values
  batchCorr_mat <- removeBatchEffect(x = mat, batch = as.factor(dataset$samples[match(colnames(mat), dataset$samples$sample_id), batch][[1]]))
  
  # Convert matrix to data.frame
  batchCorr_df <- reshape2::melt(batchCorr_mat)
  colnames(batchCorr_df) <- c("peptide_id", "sample_id", col)
  
  # Remove column with original values
  dataset_batchCorr <- dataset$peptides
  dataset_batchCorr[[col]] <- NULL
  final <- as_tibble(merge(dataset_batchCorr, batchCorr_df, all.x = T, all.y = F))
  
  return(final)
}

#' @description
#' Generates ggplot objects for a specific factor (both regular and standardized), as well as returns the data.frame with the standardized factor values
#' @param df Data.frame
#' @param x Column in data.frame to plot in the x-axis
#' @param factor Factor to standardize and plot
#' @returns List with ggplot objects and data.frame with standardized factor values
samples_scatterPlot <- function(df, x, factor) {
  # Calculate mean of factor
  m <- mean(df[[factor]])
  # Calculate standard deviation of factor
  s <- sd(df[[factor]])
  # Standardize factor
  df[[paste(factor, "std", sep = "_")]] <- (df[[factor]] - m) / s
  # Determine outliers (smaller than -3 or bigger than 3)
  df[[paste(factor, "std", "outlier", sep = "_")]] <- df[[paste(factor, "std", sep = "_")]] < -3 | df[[paste(factor, "std", sep = "_")]] > 3

  # Plot regular statistics
  lower_limit <- median(df[[factor]]) - 3 * mad(df[[factor]])
  upper_limit <- median(df[[factor]]) + 3 * mad(df[[factor]])
  p <- ggscatter(df, x = x, y = factor, color = "batch", shape = "type", xlab = "", ylab = factor, title = sprintf("Sample statistic %s", factor)) +
    labs(subtitle = sprintf("Mean(%s)=%0.5f and std(%s)=%0.5f", factor, m, factor, s)) +
    # scale_y_continuous(breaks = round(seq(min(c(lower_limit, df[[factor]])), max(c(upper_limit, df[[factor]])), by = 1))) +
    theme(axis.text.x = element_blank()) + geom_hline(yintercept = c(lower_limit, upper_limit), linetype = "dashed", color = "red")

  # Plot standardized plot
  p_standard <- ggscatter(df, x = x, y = paste(factor, "std", sep = "_"), color = "batch", shape = "type", xlab = "", ylab = sprintf("Standardized %s", factor), title = sprintf("Standardized sample statistic %s", factor)) +
    labs(subtitle = sprintf("Mean(%s)=%0.5f and std(%s)=%0.5f", factor, m, factor, s)) +
    scale_y_continuous(breaks = round(seq(min(c(-3, df[[paste(factor, "std", sep = "_")]])), max(c(3, df[[paste(factor, "std", sep = "_")]])), by = 1))) +
    theme(axis.text.x = element_blank()) + geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "red")
  return(list(p = p, p_standard = p_standard, df = df))
}

#' @description
#' Plot DIA-NN statistics in target directory
#' @param diann_stats_file Path to the file with DIA-NN statistics
#' @param metadata data.frame with columns used to indicate which samples are associated with which batches and sample types
#' @param factors Factors to plot the scatter plot for
#' @param outDir Output directory where the figures are plotted
#' @returns A directory with plots of DIA-NN statistics
qc_diann_figures <- function(diann_stats_file, metadata, factors, outDir) {
  require(readr)
  require(ggpubr)

  metadata$File.Name <- as.factor(basename(metadata$File.Name))
  diann.stats <- as.data.frame(read_tsv(file.path(diann_stats_file), comment = ""))
  diann.stats$File.Name <- basename(diann.stats$File.Name)
  diann.stats <- merge(diann.stats, metadata)

  for (factor in factors) {
    tmp <- samples_scatterPlot(df = diann.stats, x = "File.Name", factor = factor)
    png(file.path(outDir, sprintf("Figure%da_%s.png", which(factors %in% factor), factor)))
    print(tmp$p)
    graphics.off()

    png(file.path(outDir, sprintf("Figure%db_%s_standard.png", which(factors %in% factor), factor)))
    print(tmp$p_standard)
    graphics.off()

    diann.stats <- tmp$df
  }
  return(diann.stats[, c("File.Name", grep("outlier", colnames(diann.stats), value = T))])
}

#' @description
#' Convert data.frame to matrix based on specific formula using a particular column as values
#' @param df Data.frame
#' @param formula Formula to indicate which factors to use as rows and columns
#' @param col Column from the data.frame with values to put in matrix
#' @returns Matrix based on original data.frame
dfToMat <- function(df, formula = peptide_id ~ sample_id, col) {
  # Convert data.frame of intensities for peptides in samples to expression matrix
  mat <- reshape2::dcast(df, formula = formula, value.var = col)
  rownames(mat) <- mat[, 1]
  mat <- as.matrix(mat[, -1])
  return(mat)
}

#' @description
#' Run technical quality control on samples
#' @param dataset Dataset object from msdap
#' @param outDir Output directory where the figures are plotted
#' @param intensity Which intensity column to use (default is "intensity", can be used for "intensity_by_group" for noemalized data)
#' @param batch Factor used to select the sample colors in the plot (usually ms_batch)
#' @param type Factor used to select the sample shape in the plot (usually sampletype)
#' @returns A directory with the QC plots for all samples, including a table where certain statistics were compared to show whether samples are outliers or not
qc_msdap_figures <- function(dataset, outDir = NULL, intensity = "intensity", type = "sampletype", batch = "ms_batch") {
  require(ggpubr)
  require(ComplexHeatmap)
  require(circlize)
  figs <- list()

  # Set sample IDs as factors for ordered plotting
  dataset$samples$sample_id <- as.factor(dataset$samples$sample_id)
  dataset$peptides$sample_id <- as.factor(dataset$peptides$sample_id)

  # Figure 1: plot CDF of NAs in samples ------------------------------------
  mat <- dfToMat(df = dataset$peptides, col = intensity)
  df <- data.frame(sample_id = dataset$samples$sample_id, batch = dataset$samples[[batch]], type = dataset$samples[[type]])
  df$percNA <- sapply(df$sample_id, function(s) {
    tmp <- mat[, colnames(mat) %in% s]
    100 * sum(is.na(tmp)) / length(tmp)
  })
  df$percNA_33Thresh <- df$percNA > 33
  figs$Figure1_cdfNA <- ggecdf(df, x = "percNA", title = "CDF of samples with % of NAs", xlab = "% of NAs", ylab = "quantile of samples", color = "red") + labs(subtitle = "Threshold: 33% of NAs") + geom_vline(xintercept = 33, linetype = "dashed", color = "black")

  # Figure 2: standardized percentages of NAs in different samples ------------------------------------
  tmp <- try(samples_scatterPlot(df = df, x = "sample_id", factor = "percNA"), silent = T)
  if (!grepl("error", class(tmp))) {
    figs$Figure2a_percNA <- tmp$p
    figs$Figure2b_percNA_standardized <- tmp$p_standard
    df <- tmp$df
  }

  # Figure 3: standardized Kolmogorov-Smirnov distances between empirical CDF of sample precursor distribution and CDF of random normal distribution ------------------------------------
  f0 <- rowMeans(mat, na.rm = T)
  ks <- try(sapply(1:ncol(mat), function(i) {
    tmp <- ks.test(x = mat[, i], y = f0)
    return(tmp$statistic)
  }), silent=T)

  if (!grepl("error", class(ks))) {
    df$ks <- ks
    tmp <- samples_scatterPlot(df = df, x = "sample_id", factor = "ks")
    figs$`Figure3a_Kolmogorov-Smirnov_distances` <- tmp$p
    figs$`Figure3b_Kolmogorov-Smirnov_distances_standard` <- tmp$p_standard
    df <- tmp$df
  }

  # Figures 4 and 5: median and IQR intensities ------------------------------------
  for (stats in c("median", "IQR")) {
    df[[stats]] <- apply(mat, 2, get(stats), na.rm = T)
    tmp <- try(samples_scatterPlot(df = df, x = "sample_id", factor = stats), silent = T)
    if (!grepl("error", class(tmp))) {
      figName <- if (stats == "median") {
        "Figure4a_medianIntensities.png"
      } else {
        "Figure5a_iqrIntensities.png"
      }
      figs[[figName]] <- tmp$p

      figName <- if (stats == "median") {
        "Figure4b_medianIntensities_standard.png"
      } else {
        "Figure5b_iqrIntensities_standard.png"
      }
      figs[[figName]] <- tmp$p_standard

      df <- tmp$df
    }
  }

  # Figure 6: plot sample correlations ------------------------------------
  ## Retrieve sample correlations
  sample_cor <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  if (!any(apply(is.na(sample_cor), 1, all))) {
    # standardize correlations per sample
    sample_cor_std <- Reduce("rbind", lapply(1:nrow(sample_cor), function(i) {
      m <- mean(sample_cor[i, ], na.rm = T)
      s <- sd(sample_cor[i, ], na.rm = T)
      (sample_cor[i, ] - m) / s
    }))
    # calculate average standardized correlations
    sample_cor_std_avg <- apply(sample_cor_std, 2, mean, na.rm = T)
    df$avg_sample_corr <- sample_cor_std_avg[match(df$sample_id, names(sample_cor_std_avg))]

    ## Get heatmap annotations
    # All samples
    batch_vec <- df$batch[match(colnames(mat), df$sample_id)]
    b <- sort(unique(batch_vec))
    type_vec <- df$type[match(colnames(mat), df$sample_id)]
    t <- sort(unique(type_vec))
    col <- list(batch = rainbow(length(b)), type = rainbow(length(t)))
    names(col$batch) <- b
    names(col$type) <- t
    column_ha <- HeatmapAnnotation(batch = batch_vec, type = type_vec, col = col)
    ha_min <- min(sample_cor, na.rm = T)
    ha_max <- max(sample_cor, na.rm = T)
    ha_med <- median(sample_cor, na.rm = T)
    col_fun <- colorRamp2(c(ha_min, ha_med, ha_max), c("blue", "white", "darkgreen"))

    ## Plot heatmap
    figs$Figure6_sample_correlations <- Heatmap(sample_cor,
      column_title = "Sample Pearson correlation matrix",
      col = col_fun, na_col = "lightgrey", name = "Stand corr",
      clustering_distance_columns = "pearson", clustering_method_columns = "ward.D2",
      clustering_distance_rows = "pearson", clustering_method_rows = "ward.D2",
      top_annotation = column_ha, show_row_names = F, show_column_names = F
    )

    # Figure 6a: average sample correlations
    m <- mean(df$avg_sample_corr, na.rm = T)
    s <- sd(df$avg_sample_corr, na.rm = T)
    df$avg_sample_corr <- (df$avg_sample_corr - m) / s
    df$avg_sample_corr_outlier <- df$avg_sample_corr < -3 | df$avg_sample_corr > 3
    figs$Figure6a_avgSampCorr <- ggscatter(df, x = "sample_id", y = "avg_sample_corr", color = "batch", shape = "type", xlab = "", ylab = "Standardized average sample correlations", title = "Average standardized sample correlations") +
      labs(subtitle = sprintf("Mean(AvgSampCorr)=%.4e and std(AvgSampCorr)=%0.5f", m, s)) +
      scale_y_continuous(breaks = round(seq(min(c(-3, df$avg_sample_corr)), max(c(3, df$avg_sample_corr)), by = 1))) +
      theme(axis.text.x = element_blank()) + geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "red")
  }

  # Figure 7: density plot ------------------------------------
  dataset$peptides$batch <- dataset$samples[[batch]][match(dataset$peptides$sample_id, dataset$samples$sample_id)]
  dataset$peptides$sample_type <- dataset$samples[[type]][match(dataset$peptides$sample_id, dataset$samples$sample_id)]
  for (st in unique(dataset$peptides$sample_type)) {
    tib_plot <- dataset$peptides[dataset$peptides$sample_type %in% st, ]
    # plot limits
    g_fc_xlim <- abs(quantile(tib_plot[[intensity]], probs = c(.005, .995), na.rm = T))
    figs[[paste("Figure7_intensity_density_sampleType_", st, sep = "")]] <- ggplot(tib_plot, aes_string(x = intensity, labels = "sample_id", color = "batch")) +
      geom_line(stat = "density", adjust = 1, bw = "sj", na.rm = T, alpha = 0.33, linewidth = .75) +
      # scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed")) +
      xlim(g_fc_xlim) +
      labs(x = paste("log2(", intensity, ")", sep = ""), y = "Intensity density", title = sprintf("Intensity distribution by batch for sample type %s", st), colour = "Sample batch") +
      theme_bw()
  }

  # Figure 8: coefficient of variation (CoV) ------------------------------------
  # Set colors for different samples
  dataset$samples$cov_col <- paste(dataset$samples[[batch]], dataset$samples[[type]], sep = "_")
  samples_colors_long <- sample_color_coding__long_format(dataset$samples)
  samples_colors <- samples_colors_long %>%
    dplyr::select(sample_id, shortname, prop, clr) %>%
    pivot_wider(id_cols = c(sample_id, shortname), names_from = prop, values_from = clr)

  # Create coefficient of variation plot
  ## Peptides
  pep_cov <- ggplot_coefficient_of_variation(tib_input = dataset$peptides, samples = dataset$samples, samples_colors = samples_colors, column = "cov_col", intensity_col = intensity)
  figs$Figure8_CoV <- pep_cov + labs(title = "Coefficient of variation across sample types")

  if (!is.null(outDir)) {
    for (i in 1:length(figs)) {
      name <- names(figs)[i]
      p <- figs[[i]]
      png(file.path(outDir, paste(name, "png", sep = ".")))
      print(p)
      graphics.off()
    }
  }
  return(df)
}

#' @param dataset Dataset object from msdap
#' @param impute Imputation method to be used ("ppca" for probablistic PCA or "missMDA" imputePCA function)
#' @param features Create the matrix from peptide intensities or summarized protein intensities
#' @param rollup_algo If the PCA is calculated for proteins (features="proteins"), select which rollup algorithm to use
#' @param intensity Which intensity column to use (default is "intensity", can be used for "intensity_by_group" for normalized data)
#' @param quali_vars Qualitative variables used to calculate variable centres PCA plot
#' @param quanti_vars Quantitative variables used to calculate the loadings for PCA plot
#' @param calc.nPCs Number of principle components to calculate (default: 15)
#' @returns A PCA object to be plotted later, and the optimal number of PCs used to calculate the PCA (also used as the number of dimensions to include in the Scree plot)
pca_msdap_res <- function(dataset, impute = c("ppca", "missMDA"), features = c("peptides", "proteins"), rollup_algo = NULL, intensity = "intensity", quali_vars = NULL, quanti_vars = NULL, calc.nPCs = 15) {
  require(FactoMineR)
  if (features == "peptides") {
    mat <- as.data.frame(dfToMat(dataset$peptides, col = intensity, formula = sample_id ~ peptide_id))
  } else if (features == "proteins") {
    if (!is.null(rollup_algo)) {
      tmp <- dataset$peptides[, c("peptide_id", "protein_id", "sample_id", intensity)]
      colnames(tmp) <- c("peptide_id", "protein_id", "sample_id", "intensity")
      mat <- t(rollup_pep2prot(tmp, rollup_algorithm = rollup_algo, return_as_matrix = TRUE))
    } else {
      stop("No rollup algorithm selected")
    }
  } else {
    stop("No features selected. Select either peptides or proteins")
  }
  mat_scaled <- base::scale(as.matrix(mat), center = TRUE, scale = FALSE)
  
  if (length(quanti_vars) > 0) {
    mat <- merge(mat_scaled, dataset$samples[, c("sample_id", quanti_vars)], by.x = "row.names", by.y = "sample_id")
    rownames(mat) <- mat$Row.names
    mat_scaled <- mat[, -1]
  }
  
  # Remove any quantitative variables that are completely 0
  empty_vars <- names(which(colSums(mat_scaled[, quanti_vars], na.rm = T) == 0))
  quanti_vars <- setdiff(quanti_vars, empty_vars)
  mat_scaled <- mat_scaled[, !(colnames(mat_scaled) %in% empty_vars)]
  
  cat("Impute data used for PCA calculation\n")
  start <- Sys.time()
  dfcompleted <- if (impute == "ppca") {
    pcaMethods::pca(mat_scaled, method = "ppca", nPcs = calc.nPCs, seed = 123, maxIterations = 2000)@completeObs
  } else if (impute == "missMDA") {
    missMDA::imputePCA(mat_scaled, ncp = calc.nPCs, quanti.sup = if (length(quanti_vars) > 0) {
      which(colnames(mat_scaled) %in% quanti_vars)
    } else {
      NULL
    })$completeObs
  } else {
    stop("Imputation method selected is invalid")
  }
  end <- Sys.time()
  print(end - start)
  
  # Estimate optimal number of PCs while excluding quantifiable variables (should not contribute to optimal number of PCs)
  cat("Estimate optimal number of PCs\n")
  start <- Sys.time()
  estim_pcs <- estim_ncp(dfcompleted[, !(colnames(dfcompleted) %in% quanti_vars)], ncp.max = calc.nPCs, scale = T)
  ncp <- if(estim_pcs$ncp == 0){calc.nPCs-2}else{estim_pcs$ncp}
  end <- Sys.time()
  print(end - start)
  
  if (length(quali_vars) > 0) {
    mat <- merge(dfcompleted, dataset$samples[, c("sample_id", quali_vars)], by.x = "row.names", by.y = "sample_id")
    rownames(mat) <- mat$Row.names
    dfcompleted <- mat[, -1]
  }
  
  # Calculate PCA
  cat("Calculating PCA\n")
  start <- Sys.time()
  res.pca <- PCA(dfcompleted,
                 scale.unit = TRUE, ncp = ncp, ind.sup = NULL,
                 quali.sup = if (length(quali_vars) > 0) {
                   which(colnames(dfcompleted) %in% quali_vars)
                 } else {
                   NULL
                 },
                 quanti.sup = if (length(quanti_vars) > 0) {
                   which(colnames(dfcompleted) %in% quanti_vars)
                 } else {
                   NULL
                 },
                 graph = FALSE
  )
  end <- Sys.time()
  print(end - start)
  
  return(list(res.pca = res.pca, calc.nPCs = ncp))
}

#' Get optimal number of PCs to use
#' @param res.pca PCA object
#' @returns The optimal number of PCs
pca_optimal <- function(res.pca) {
  require(factoextra)
  # require(bigpca)
  # Two main methods to determine the optimal number of PCs to plot
  ## Elbow method: in the scree plot, find when a substantial drop occurs in the variance explained / eigenvalues
  ## quick.elbow function in bigpca R package can calculate it numerically
  eig <- as.data.frame(get_eig(res.pca))
  eig_opt <- quick.elbow(eig$eigenvalue, low = 0.08, max.pc = 0.9)
  var_opt <- quick.elbow(eig$variance.percent, low = 0.08, max.pc = 0.9)
  
  ## Kaiser's Rule: According to Kaiser’s rule, all the components with eigenvalues greater than 1 should be kept for statistical inference
  kaiser_opt <- nrow(eig[eig$eigenvalue > 1, ])
  
  # nb <- missMDA::estim_ncpPCA(mat_scaled, ncp.max=calc.nPCs) ## Time consuming, nb = 2)
  # calc.nPCs <- nb$ncp
  plot.nPCs <- min(kaiser_opt, min(eig_opt, var_opt))
  
  return(plot.nPCs)
}

#' @param res.pca PCA object
#' @param outDir Output directory where the figures are plotted
#' @param type Column use to separate the samples into groups
#' @param calc.nPCs Number of principle components to calculate (default: 10)
#' @param plot.nPCs Number of principle components to plot (default: 3)
#' @returns A directory with the PCA plots with all samples, separated by batch and type. Potentially includes loadings of technical and/or clinical variables
pca_msdap_figures <- function(res.pca, outDir, type = "ms_batch", calc.nPCs = 10, plot.nPCs = 3) {
  setwd(outDir)
  
  df.eig <- data.frame(dim = factor(1:nrow(get_eig(res.pca))), get_eig(res.pca))[1:calc.nPCs, ]
  
  png("Figure9a_Scree_eigenvalues.png")
  print(ggpubr::ggpar(ggplot(df.eig, aes(dim, eigenvalue, group = 1)) +
                        geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
                        geom_line(color = "black", linetype = "solid") +
                        geom_point(shape = 19, color = "black") +
                        geom_text(
                          label = round(df.eig$eigenvalue, 1), vjust = -0.4,
                          hjust = 0
                        ) +
                        labs(title = "Scree plot eigenvalues", x = "Dimensions", y = "eiganvalue")), ggtheme = theme_minimal())
  graphics.off()
  
  png("Figure9a_Scree_var_explained.png")
  print(ggpubr::ggpar(ggplot(df.eig, aes(dim, variance.percent, group = 1)) +
                        geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
                        geom_line(color = "black", linetype = "solid") +
                        geom_point(shape = 19, color = "black") +
                        geom_text(
                          label = paste0(round(df.eig$variance.percent, 1), "%"), vjust = -0.4,
                          hjust = 0
                        ) +
                        labs(title = "Scree plot variances explained", x = "Dimensions", y = "variance percentage")), ggtheme = theme_minimal())
  graphics.off()
  
  png("Figure9a_Scree_cumulative_var_explained.png")
  print(ggpubr::ggpar(ggplot(df.eig, aes(dim, cumulative.variance.percent, group = 1)) +
                        geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
                        geom_line(color = "black", linetype = "solid") +
                        geom_point(shape = 19, color = "black") +
                        geom_text(
                          label = paste0(round(df.eig$cumulative.variance.percent, 1), "%"), vjust = -0.4,
                          hjust = 0
                        ) +
                        labs(title = "Scree plot cumulative variance explained", x = "Dimensions", y = "cumulative variance percentage")), ggtheme = theme_minimal())
  graphics.off()
  
  # Plot PCA
  for (i in 1:(plot.nPCs - 1)) {
    for (j in (i + 1):plot.nPCs) {
      # Plot ellipses around qualitative variables
      p <- try(plotellipses(res.pca, axes = c(i, j), level = 0.99, means = FALSE), silent = T)
      if (any(!grepl("error", class(p)))) {
        png(paste("Figure9b_PCA_ellipses_qualiSupp_PC", i, "_vs_PC", j, ".png", sep = ""))
        print(p)
        graphics.off()
      }
      # Plot
      png(paste("Figure9c_PCA_qualitativeCenters_PC", i, "_vs_PC", j, ".png", sep = ""))
      print(plot.PCA(res.pca, choix = "ind", axes = c(i, j), invisible = "ind", title = paste("Plot centers of quantitative variables PC", i, " vs PC", j, sep = "")))
      graphics.off()
      # Plot biplot of samples based on batches with supplementary quantitative variables
      png(paste("Figure9d_PCA_factoexra_biplot_PC", i, "_vs_PC", j, ".png", sep = ""))
      print(fviz_pca_biplot(X = res.pca, axes = c(i, j), geom.ind = c("point"), habillage = type, addEllipses = TRUE, ellipse.level = 0.99, title = paste("Plot individual samples PC", i, " vs PC", j, sep = ""), invisible = "var", col.quanti.sup = "blue"))
      graphics.off()
    }
  }
}

#' @description
#' fast CoV computation for all rows in a matrix of natural-log transformed values
#' @param matrix_with_naturallog_intensities a matrix with values transformed into natural log
#' @importFrom matrixStats rowSds
coefficient_of_variation_vectorized <- function(matrix_with_naturallog_intensities) {
  sqrt(expm1(matrixStats::rowSds(matrix_with_naturallog_intensities, na.rm = T)^2))
}

#' @description
#' wrapper function around tidyverse's gather to convert a matrix to long-format tibble
#' NA values are removed from output
#' @param mat a matrix
#' @param value_name target column name for the values
#' @param column_name target column name for the row names
#' @param row_name target column name for the column names
#' @importFrom tidyr gather
matrix_to_long <- function(mat, value_name = "value", column_name = "sample", row_name = "sequence") {
  x <- as_tibble(mat)
  x$rownms <- rownames(mat)
  res <- x %>% gather(colname, value, -rownms, na.rm = T)
  names(res)[1] <- row_name
  names(res)[2] <- column_name
  names(res)[3] <- value_name
  return(res)
}

#' @description
#' Plot coefficient of variations (CoVs) in different groups
#' @param tib_input peptides from dataset
#' @param samples samples from dataset
#' @param samples_colors colors for metadata
#' @param column column used to plot CoV by
#' @param intensity_col Column with intensity values
#'
#' @importFrom ggpubr theme_pubr
#' @importFrom matrixStats rowSums2
ggplot_coefficient_of_variation <- function(tib_input, samples, samples_colors, column, intensity_col = "intensity") {
  ## for CoV computation, we need natural log while intensities are log2; log_b(x) = log_d(x) / log_d(b)  @  https://www.purplemath.com/modules/logrules5.htm
  # toy example; given y = log2(x=100), we need z = log10(x);
  # x = 100
  # y = log2(x)
  # z = y / log2(10)
  # x;y;z

  # !! here we use the by-group filtering, and normalization, as configured by the user
  # this is important, because minpep, topN and different normalization make it different from leave-one-out CoV (there, we must use 'mode' norm for speed as we have to normalize the dataset as often as there are samples)
  tib_input <- tib_input[, c("peptide_id", "sample_id", intensity_col)]
  colnames(tib_input) <- c("peptide_id", "sample_id", "intensity")
  tibw_abundance_naturallog <- tib_input %>%
    filter(!is.na(intensity)) %>% # this drops all peptides not passing the filter in any sample/group, which makes downstream wide format intensity matrix much smaller
    mutate(intensity = intensity / log2(exp(1))) %>%
    pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity)

  ## compute stats per group
  groups <- unique(samples[[column]])
  mat_cov <- matrix(NA, nrow = nrow(tibw_abundance_naturallog), ncol = length(groups), dimnames = list(tibw_abundance_naturallog$peptide_id, groups))
  dropcols <- NULL
  for (grp in groups) {
    # samples for current group
    sid <- intersect(samples$sample_id[samples[[column]] %in% grp], colnames(tibw_abundance_naturallog))
    # skip if <n samples in group
    if (length(sid) < 2) {
      dropcols <- c(dropcols, grp)
      append_log(sprintf("No CoV computation for sample group '%s', require at least 2 replicates", grp), type = "info")
      next
    }

    # fast CoV computation
    m <- as.matrix(tibw_abundance_naturallog %>% dplyr::select(!!sid))
    rows_fail <- matrixStats::rowSums2(!is.na(m)) < 3
    # less than 50 peptides have a value, not a meaningful set of data.points for CoV analysis
    if (sum(!rows_fail) < 50) {
      dropcols <- c(dropcols, grp)
      append_log(sprintf("No CoV computation for sample group '%s', requires at least 50 peptides/proteins to have a value", grp), type = "info")
      next
    }
    # m[rows_fail, ] = NA # remove rows with less than 3 values (can technically calculate on 2 values, but we chose to require at least 3)
    mat_cov[!rows_fail, grp] <- coefficient_of_variation_vectorized(m[!rows_fail, ])
  }

  # remove sample groups that have too few replicates
  mat_cov <- mat_cov[, !(colnames(mat_cov) %in% dropcols), drop = F]
  if (ncol(mat_cov) == 0) {
    append_log("No data available for CoV computation, skipping plots", type = "info")
    return(list())
  }

  # cov as percentages in all downstream analyses
  mat_cov <- mat_cov * 100
  # debug/QC: bp=boxplot(mat_cov, outline=F, las=2); boxplot_add_text(bp, cex=.5)

  ## summary stats: boxplot. we plot these as text labels onto the ggplot downstream
  bp <- boxplot(mat_cov, plot = F)
  cov_data_summ <- data.frame(bp$stats)
  colnames(cov_data_summ) <- bp$names
  # We split the text labels for the CoV boxplot (all minus last  vs  last) so we can apply different vertical justification
  cov_data_summ_long_14 <- cov_data_summ[c(1, 4), , drop = F] %>%
    tidyr::gather(group, summ) %>%
    dplyr::rename(cov = summ)
  cov_data_summ_long_23 <- cov_data_summ[c(2, 3), , drop = F] %>%
    tidyr::gather(group, summ) %>%
    dplyr::rename(cov = summ)
  cov_data_summ_long_5 <- cov_data_summ[5, , drop = F] %>%
    tidyr::gather(group, summ) %>%
    dplyr::rename(cov = summ)

  ## summary stats: number of peptides used in CoV computation per group
  tib_cov_group_n <- tibble::enframe(colSums(!is.na(mat_cov)), name = "group", value = "n")
  # tib_cov_group_n$cov_median = as.numeric(cov_data_summ[3, match(colnames(cov_data_summ), tib_cov_group_n$group)]) # coordinates for plotting at boxplot median line

  # adjust text size to the number of groups
  n_groups <- ncol(mat_cov)
  txt_size <- 4 - 2 * min(10, n_groups) / 20

  # prepare plot data in long format
  tib_plot <- matrix_to_long(mat_cov, value_name = "cov", column_name = "group", row_name = "peptide_id") %>%
    mutate(group = factor(group, levels = colnames(mat_cov)))

  ## boxplot
  p_boxplot <- ggplot(tib_plot, aes(group, cov, fill = group)) +
    geom_boxplot(outlier.shape = NA, na.rm = T) +
    geom_text(data = cov_data_summ_long_14, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = -.5, nudge_x = .05, size = txt_size) +
    geom_text(data = cov_data_summ_long_23, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0.5, vjust = -.5, size = txt_size) +
    geom_text(data = cov_data_summ_long_5, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = 1, nudge_x = .05, size = txt_size) +
    geom_text(data = tib_cov_group_n, aes(x = group, y = -4, label = n), hjust = 0.5, vjust = 0.5, size = txt_size) +
    # geom_text(data = cov_data_summ_long_2, aes(x = group, y = cov, label = round(cov, digits = 1)), hjust = 0, vjust = .5, nudge_x = .05, size = txt_size) +
    # geom_text(data = tib_cov_group_n, aes(x = group, y = cov_median, label = paste0("n=",n)), hjust = 1, vjust = -.5, nudge_x = -.05, size = txt_size) +
    scale_color_manual(values = array(unique(samples_colors[[column]]), dimnames = list(unique(samples[[column]]))), aesthetics = c("fill")) +
    coord_cartesian(ylim = c(-5, max(cov_data_summ)), clip = "off") +
    scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(cov_data_summ, na.rm = T) / 10) * 10, by = 10)) +
    labs(x = "", y = "Coefficient of Variation (%)") +
    ggpubr::theme_pubr() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = ifelse(n_groups < 8, 11, 8)))
  return(p_boxplot)
}
