library(jsonlite)
library(data.table)

#' @title Stop with a clear message if IDs (row names) are duplicated
#' @description
#' read.table(row.names = 1) aborts with the cryptic
#' "duplicate 'row.names' are not allowed" and does not tell you WHICH file or
#' WHICH ids are the problem. This helper is called on the intended row-name
#' vector BEFORE it is used, so we can raise a message that actually names the
#' source file, the column, and the offending ids.
#' @param ids Character vector of intended row names (e.g. feature or sample ids).
#' @param source_file Path of the file the ids came from (for the message).
#' @param id_label Human-readable name of the id column (e.g. "Run" or "feature").
#' @return Invisibly TRUE if all ids are unique; otherwise stops execution.
check_unique_ids <- function(ids, source_file, id_label) {
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    preview <- paste(utils::head(dup, 10), collapse = ", ")
    if (length(dup) > 10) preview <- paste0(preview, ", ...")
    stop(sprintf(
      paste0("Duplicate %s ids found in '%s' (%d duplicated value(s)): %s\n",
             "  Row/sample ids must be unique. For DIA-NN precursors this",
             " usually means the channel was not added to the label at",
             " every level (precursor/peptide/protein). Add the channel to",
             " the name, apply a filter (e.g. Channel == L), or de-duplicate",
             " before this step."),
      id_label, source_file, length(dup), preview),
      call. = FALSE)
  }
  invisible(TRUE)
}

#' @title Read a TSV file the same way pandas does in the Python modules
#' @description
#' Reads a tab-separated file with the quoting rules used by
#' \code{pandas.read_csv(sep = "\t")}: only \code{"} is a quote character
#' (\code{""} inside a quoted field is a literal \code{"}), there is no comment
#' character, so \code{'} and \code{#} are read as plain text. Unlike
#' \code{read.table}, ids containing these characters are not truncated or
#' merged across rows.
#' @param path Path to the file.
#' @param row_names If TRUE, the first column is used as row names
#'   (like \code{index_col = 0} in pandas) after checking it for duplicates.
#' @param id_label Name of the id column used in the duplicate error message.
#' @param sep Field separator, "\t" by default ("auto" lets fread detect it).
#' @return data.frame; with row_names = TRUE the first column becomes row names.
read_tsv_pd <- function(path, row_names = TRUE, id_label = "feature", sep = "\t") {
  df <- fread(path,
              sep = sep,
              quote = "\"",
              header = TRUE,
              skip = 0,
              na.strings = c("", "NA", "NaN"),
              integer64 = "double",
              keepLeadingZeros = TRUE,
              check.names = FALSE,
              data.table = FALSE)
  # fread strips the outer quotes of a quoted field but keeps the escaped ""
  # inside it, so turn them back into " (as pandas does)
  unescape <- function(x) gsub("\"\"", "\"", x, fixed = TRUE)
  names(df) <- unescape(names(df))
  chr_cols <- vapply(df, is.character, logical(1))
  df[chr_cols] <- lapply(df[chr_cols], unescape)
  if (row_names) {
    ids <- as.character(df[[1]])
    check_unique_ids(ids, path, id_label)
    df[[1]] <- NULL
    rownames(df) <- ids
  }
  df
}

#' @title Write a TSV file the same way pandas does in the Python modules
#' @description
#' Writes row names as the first column (like \code{DataFrame.to_csv(sep = "\t")})
#' and quotes only fields containing a tab, \code{"} or a newline, doubling
#' inner \code{"} (pandas' \code{QUOTE_MINIMAL}). Missing values are written as
#' empty fields.
#' @param df data.frame or matrix to write.
#' @param path Output file path.
#' @return Invisibly NULL.
write_tsv_pd <- function(df, path) {
  fwrite(as.data.frame(df),
         path,
         sep = "\t",
         quote = "auto",
         qmethod = "double",
         row.names = TRUE,
         na = "")
  invisible(NULL)
}

#' @title Run R Code with Logging (Append to Existing Log)
#' @description
#' Appends all console output, messages, and warnings from an R expression to an existing log file.
#' Errors are not caught — they stop the script, ensuring pipeline failure detection.
#'
#' @param expr The R expression to evaluate (e.g. your main processing call).
#' @param log_file Path to an existing log file where new output should be appended.
#'
#' @return The result of evaluating \code{expr}, or stops if an error occurs.
run_with_logging <- function(expr, log_path) {
  
  con <- file(log_path, open = "at")  # append mode
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  
  ## Add a time-stamped header 
  cat(sprintf("\n\n--- Module in R -- START [%s] ---\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
      file = con)
  
  # Redirect both standard output and messages to log file
  sink(con, type = "output")
  sink(con, type = "message")
  
  result <- tryCatch(
    withCallingHandlers(
      expr(),
      warning = function(w) {
        message("###WARNING### : ", conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        message("MESSAGE: ", conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      ## log the real error, then re-raise so the module halts with it
      ## (returning NULL here masks the cause: downstream code crashes on the
      ##  NULL result with an unrelated message instead of this one)
      message("###ERROR###: ", conditionMessage(e))
      stop(e)
    }
  )
  
  ## Mark end of block
  cat(sprintf("\n--- Module in R -- END [%s] ---\n", 
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
      file = con)
  
  return(result)
}

#' @title Read JSON file and convert empty lists to NULL recursively
#'
#' @description This function reads a JSON file into an R list using 
#' \code{jsonlite::fromJSON}.It then recursively converts any empty list 
#' elements (i.e., \code{list()}) into \code{NULL}, which helps with consistent
#' handling of missing or empty values.
#'
#' @param parameter_file Character string. Path to the JSON parameter file.
#' 
#' @return The parameter file where empty lists are replaced with \code{NULL}.
read_JSON <- function(parameter_file){
  
  ## read in parameter file
  parameter_list <- jsonlite::fromJSON(parameter_file)
  
  # convert empty list elements to null
  emptyListToNull <- function(x) {
    if (is.list(x) && length(x) == 0) {
      return(NULL)
    } else if (is.list(x)) {
      return(lapply(x, emptyListToNull))
    } else {
      return(x)
    }
  }
  parameter_list <- emptyListToNull(parameter_list)
  
  return(parameter_list)
} 

#' @title 
#' Create a command-line argument parser for the R modules
#'
#' @desciption
#' Constructs an \code{optparse::OptionParser} object that defines the
#' command-line interface for the R scripts. The parser specifies required and 
#' optional arguments such as the input data directory, parameter file, output 
#' directory, metadata file, IDs file, step counter, and log file.
#'
#' @return An \code{OptionParser} object, ready to be passed to
#'   \code{optparse::parse_args()}.
#'
#' @details
#' The following options are supported:
#' \itemize{
#'   \item \code{-d, --data}: Path to the input data file, a .tsv file (required).
#'   \item \code{-p, --params}: Path to a JSON file with normalization parameters (required).
#'   \item \code{-o, --out}: Path to the output directory (required).
#'   \item \code{-m, --md}: Path to a metadata file (optional).
#'   \item \code{-i, --ids}: Path to an IDs file (optional).
#'   \item \code{-c, --counter}: Pipeline step counter as an integer (required).
#'   \item \code{-l, --log}: Path to a log file (required).
#' }
make_parser <- function() {
  option_list <- list(
    make_option(c("-d","--data"),     type="character", help="Path to data file"),
    make_option(c("-p","--params"),   type="character", help="Path to params JSON"),
    make_option(c("-o","--out"),      type="character", help="Output directory"),
    make_option(c("-m", "--md"), type="character", default = NULL, help="Path to metadata file (optional)"),
    make_option(c("-i","--ids"),      type="character", default = NULL, help="Path to IDs file (optional)"),
    make_option(c("-c","--counter"),  type="character",   help="Pipeline step counter"),
    make_option(c("-l","--log"),      type="character", help="Log file path")
  )
  OptionParser(option_list = option_list,
               description = "Normalize DIA-NN output in R")
}

#' @title
#' Validate parsed command-line options for the R modules 
#'
#' @description
#' Performs validation checks on the options returned by
#' \code{optparse::parse_args()}. Ensures that required arguments are present,
#' required files or directories exist, and that output/log directories are
#' valid. Stops execution with an error message if validation fails.
#'
#' @param opt A list of options returned by \code{optparse::parse_args()}.
#'
#' @return Invisibly returns \code{NULL}. Its point is that it validates the 
#' input arguments and stops with errors when any of them is invalid.
#'
#' @details
#' Validation rules:
#' \itemize{
#'   \item \code{--data}: Must be provided and must be an existing file. 
#'   \item \code{--params}: Must be provided and must be an existing file.
#'   \item \code{--out}: Must be provided; will be created if it does not exist.
#'   \item \code{--log}: Must be provided; parent directory must exist.
#'   \item \code{--counter}: Must be provided and numeric.
#'   \item \code{--metadata}: Optional; if provided, must exist as a file.
#'   \item \code{--ids}: Optional; if provided, must exist as a file.
#' }
validate_opts <- function(opt) {
  # Check required arguments
  if (is.null(opt$data))
    stop("Missing required argument: --data", call. = FALSE)
  if (is.null(opt$params))
    stop("Missing required argument: --params", call. = FALSE)
  if (is.null(opt$out))
    stop("Missing required argument: --out", call. = FALSE)
  if (is.null(opt$log))
    stop("Missing required argument: --log", call. = FALSE)
  if (is.null(opt$counter))
    stop("Missing required argument: --counter", call. = FALSE)
  
  # Validate data file
  if (!file.exists(opt$data))
    stop(sprintf("Data file not found: %s", opt$data), call. = FALSE)
  
  # Validate params JSON file
  if (!file.exists(opt$params))
    stop(sprintf("Parameter file not found: %s", opt$params), call. = FALSE)
  
  # Validate output directory (create if missing)
  if (!dir.exists(opt$out)) {
    dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Validate log file parent directory
  log_dir <- dirname(opt$log)
  if (!dir.exists(log_dir))
    stop(sprintf("Log directory not found: %s", log_dir), call. = FALSE)
  
  # Validate counter
  if (suppressWarnings(is.na(as.numeric(opt$counter))))
    stop("--counter must be numeric", call. = FALSE)
  
  # Optional metadata file
  if (!is.null(opt$metadata) && !file.exists(opt$metadata))
    stop(sprintf("Metadata file not found: %s", opt$metadata), call. = FALSE)
  
  # Optional IDs file
  if (!is.null(opt$ids) && !file.exists(opt$ids))
    stop(sprintf("IDs file not found: %s", opt$ids), call. = FALSE)
  
  invisible(NULL)
}