library(jsonlite)

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
    
    # Evaluate the expression (e.g. main())
    result <- expr()
    
    
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
    

