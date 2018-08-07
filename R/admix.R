# admix.R
#' Constructor for admix object
#'
#' @return an admix object which is a list with 6 elements:
#'  K: number of clusters estimated by ADMIXTURE
#'
#'  nsamples: number of samples used
#'
#'  nmarkers: number of markers used
#'
#'  Q_df: a data.frame of cluster membership probabilities
#'
#'  P_df: a data.frame of estimated marker frequencies in each inferred population
#'
#'  Q_bias_df: a data.frame of bias estimates in cluster membership probabilities
#'
#'  Q_se_df: a data.frame of standard error estimates in cluster membership probabilities
#'
#'  fst_matrix: a symmetric matrix of F_st values between inferred clusters
#'
#'  log_info: a data.frame containing the K, CVerror and logLik of the last model.
#'
#' @export
admix <- function() {
  structure(list(K = NULL, nsamples = NULL, nmarkers = NULL,
                 Q_df =  NULL, P_df = NULL, Q_bias_df = NULL,
                 Q_se_df = NULL, fst_matrix = NULL, log_info = NULL),
            class = "admix")
}

#' Contructor for admixList
#'
#' @description Collect many  \code{\link{admix}} objects
#' @param ... a list of \code{\link{admix}} objects
#' @return an admixList object
#' @importFrom methods is
#' @export
admixList <- function(...) {
  listOut <- list(...)

  if (length(listOut) == 1L && is.list(listOut[[1L]])) {
    listOut <- listOut[[1L]]
    if (inherits(listOut, "admix")) {
      return(listOut)
    }

  }
  if (length(listOut) == 0L) {
    stop("an admixList must contain at least one admix object")
  } else {
    if (!all(sapply(listOut, is, 'admix')))
      stop("all elements in '...' must be admix objects")
    structure(listOut, class = "admixList")
  }
}

#' Example admixture runs
#' @export
exampleAdmixture <- function() {
  path <- system.file("extdata/hapmap3_files", package = "starmie")
  logs <- list.files(path, pattern = "*.out$", full.names = TRUE)
  qfin <- list.files(path, pattern = "*.Q$", full.names = TRUE)
  pfin <- list.files(path, pattern = "*.P$", full.names = TRUE)
  admixList(mapply(loadAdmixture, qfin, pfin, logs,
                   SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

#' @export
print.admix <- function(x, ...) {
  cat(paste("admix object containing run information for k =", x$K, "\n"))
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
print.admixList <- function(x, ...) {
  n <- length(x)
  cat(paste("admixList object containing", n, " ADMIXTURE runs.\n"))
  cat(paste("Number of Ks by number of runs:"))
  Ks <- table(unlist(lapply(x, getK)))
  print(Ks)
  cat(paste("Model fit information by K:\n"))
  log_all <- combineLogs(x)
  print(log_all)
}

combineLogs <- function(x) {
  stopifnot(inherits(x, "admixList"))
  do.call("rbind", lapply(x, function(y) y$log_info))
}

#' @export
`[[.admixList` <- function(x, ...) {

  y <- NextMethod("[[")
  y
}

#' @export
`[.admixList` <- function(x, ...) {
  y <- NextMethod("[")
  if (length(y) == 1L) {
    return(y)
  } else {
    class(y) <- oldClass(x)
    y
  }
}
