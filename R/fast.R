# fast.R
#' Constructor for fast object
#'
#' @return a fast object which is a list with 5 elements:
#'  K: number of clusters estimated by fastSTRUCTURE
#'
#'  nsamples: number of samples used
#'
#'  nmarkers: number of markers used
#'
#'  Q_df: a data.frame of cluster membership probabilities
#'
#'  P_df: a data.frame of estimated marker frequencies in each inferred population
#'
#'  log_info: a data.frame containing the K, CVerror and logLik of the last model.
#'
#' @export
fast <- function() {
  structure(list(K = NULL, nsamples = NULL, nmarkers = NULL,
                 Q_df =  NULL, P_df = NULL, log_info = NULL),
            class = "fast")
}

#' Contructor for fastList
#'
#' @description Collect many  \code{\link{fast}} objects
#' @param ... a list of \code{\link{fast}} objects
#' @return an fastList object
#' @importFrom methods is
#' @export
fastList <- function(...) {
  listOut <- list(...)

  if (length(listOut) == 1L && is.list(listOut[[1L]])) {
    listOut <- listOut[[1L]]
    if (inherits(listOut, "fast")) {
      return(listOut)
    }

  }
  if (length(listOut) == 0L) {
    stop("an fastList must contain at least one fast object")
  } else {
    if (!all(sapply(listOut, is, 'fast')))
      stop("all elements in '...' must be fast objects")
    structure(listOut, class = "fastList")
  }
}

#' Example fastStructure runs
#' @export
exampleFastStruct <- function() {
  path <- system.file("extdata/hapmap3_files", package = "starmie")
  logs <- list.files(path, pattern = "*.log$", full.names = TRUE)
  qfin <- list.files(path, pattern = "*.meanQ$", full.names = TRUE)
  pfin <- list.files(path, pattern = "*.meanP$", full.names = TRUE)
  fastList(mapply(loadFastStructure, qfin, pfin, logs,
                   SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

#' @export
print.fast <- function(x, ...) {
  cat(paste("fast object containing run information for k =", x$K, "\n"))
  cat("Model parameters:\n")
  cat(paste("  ", "No. samples:", x$nsamples, "\n"))
  cat(paste("  ", "No. markers:", x$nmarkers, "\n"))
  logfile <- !is.null(x$log_info)
  cat(paste("Model diagnostics available:", logfile, "\n"))
  if (logfile) {
    cat("Model fit statistics:\n")
    print(x$log_info)
  }

}

#' @export
print.fastList <- function(x, ...) {
  n <- length(x)
  cat(paste("fastList object containing", n, " fastSTRUCTURE runs.\n"))
  cat(paste("Number of Ks by number of runs:"))
  Ks <- table(unlist(lapply(x, getK)))
  print(Ks)
  cat(paste("Model fit information by K:\n"))
  log_all <- combineFastLogs(x)
  print(log_all)
}

combineFastLogs <- function(x) {
  stopifnot(inherits(x, "fastList"))
  do.call("rbind", lapply(x, function(y) y$log_info))
}

#' @export
`[[.fastList` <- function(x, ...) {

  y <- NextMethod("[[")
  y
}

#' @export
`[.fastList` <- function(x, ...) {
  y <- NextMethod("[")
  if (length(y) == 1L) {
    return(y)
  } else {
    class(y) <- oldClass(x)
    y
  }
}
