# starmie.R
# Place documentation for test datasets here.

#' Constructor for struct object
#' @description struct object for storing structure run information
#' @details The \code{\link{struct}} object is a list with 11 elements:
#'
#'  K: number of clusters estimated by Structure
#'
#'  run_params: the run parameters given to the Structure program
#'
#'  mem_df: assigned cluster membership proportions
#'
#'  allele_freqs: Estimated net nucleotide distances between cluster
#'
#'  avg_dist_df: Average distance between individuals
#'
#'  fit_stats_df: Model fit statistics
#'
#'  fst_df: Fst values
#'
#'  ancest_df: Inferred ancestry of individuals
#'
#'  clust_allele_list: Cluster allele frequencies
#'
#'  burn_df: Burn in MCMC iteration output
#'
#'  nonburn_df: Main MCMC iteration output
#' @return a \code{\link{struct}} object
#' @export
#' @seealso \code{\link{loadStructure}} for reading in STRUCTURE out_f files.
#'
#' \code{\link{structList}} for manipulating multiple struct objects
struct <- function() {
  structure(list(K = NULL, run_params = NULL, mem_df = NULL,
                 allele_freqs = NULL, avg_dist_df = NULL,
                 fst_df=NULL, fit_stats_df = NULL,  ancest_df = NULL,
                 clust_allele_list = NULL,
                 burn_df=NULL, nonburn_df=NULL),
            class = "struct")
}

#' Constructor for a structList object
#' @description  the structList class is a conatainer for storing a collection
#' of struct objects.
#' @param ... a list of a \code{\link{struct}} objects
#' @importFrom methods is
#' @export
structList <- function(...) {

  listOut <- list(...)

  if (length(listOut) == 1L && is.list(listOut[[1L]])) {
    listOut <- listOut[[1L]]
    if (inherits(listOut, "struct")) {
      return(listOut)
    }

  }
  if (length(listOut) == 0L) {
    stop("a structList must contain at least one struct object")
  } else {
    if (!all(sapply(listOut, is, 'struct')))
      stop("all elements in '...' must be struct objects")
    structure(listOut, class = "structList")
  }
}

#' Accessor methods for struct objects
#' @description getD Return the number of free parameters in STRUCTURE model
#' @param structure_obj a \code{\link{struct}} object
#' @export
getD <- function(structure_obj) {

  k <- getK(structure_obj)
  n <- as.integer(structure_obj$run_params[1,2])
  l <- as.integer(structure_obj$run_params[2,2])
  J <-sum(unlist(lapply(structure_obj$clust_allele_list, function(x) x$AlleleNumber)))

  return(list(d = k*(J -l), n = n))

}

#' @describeIn getD Return the estimated log posterior probability (L_k) from a \code{\link{struct}} object
#' @export
getPosterior <- function(structure_obj){
  structure_obj$fit_stats_df[structure_obj$fit_stats_df$Statistic=="Estimated Ln Prob of Data",][2]
}

#' @describeIn getD Return the estimated mean and variance of estimated log-likelihood from a \code{\link{struct}} object
#' @export
getFitStats <- function(structure_obj) {
  structure_obj$fit_stats_df[c(2,3), "Value"]
}


#' @describeIn getD Return non-burn in MCMC iterations.
#' @export
getMCMC <- function(structure_obj) {
  if(is.null(structure_obj$nonburn_df))
    stop("Log file required to get MCMC diagnsotics")

  mcmc_df <- data.frame(K = getK(structure_obj),
                        structure_obj$nonburn_df[, c("Rep#:", "Alpha", "Ln Like")])
  colnames(mcmc_df) <- c("K", "Iteration", "Alpha", "LogL")
  mcmc_df
}

#' Example structure objects
#' @param example_type a character string either "multiple_runs", "clumpp" or
#' "mcmc_diagnostics" or "barplot"
#' @description load structure objects for different starmie functions
#' @export
exampleStructure <- function(example_type) {
  structure_files <- system.file("extdata/microsat_testfiles",
                                 package="starmie")

  if (example_type == "clumpp") {
    k3_out <- list.files(structure_files,
                         pattern = ".*K3.*out_f",
                         full.names = TRUE)

    lapply(k3_out, loadStructure)
  } else if (example_type == "mcmc_diagnostics") {
    k10_out <- list.files(structure_files,
                         pattern = ".*K10.*out_f",
                         full.names = TRUE)
    k10_log <- list.files(structure_files,
                         pattern = ".*K10.*log",
                         full.names = TRUE)
    structList(mapply(loadStructure, k10_out, k10_log, SIMPLIFY = FALSE))
  } else if (example_type == "barplot") {
    k6_out <- list.files(structure_files,
                         pattern = "^locprior_K6.out_f",
                         full.names = TRUE)
    loadStructure(k6_out)

  } else if (example_type == "multiple_runs") {
    runs <-  list.files(structure_files,
                        pattern = "out_f$",
                        full.names = TRUE)
    sequential <- runs[grepl("^run2|^locprior", basename(runs))]
    structList(lapply(sequential, loadStructure))
  }

}

#' @export
print.struct <- function(x, ...) {
  cat(paste("struct object containing run information for k =", x$K, "\n"))
  cat("Model run parameters:\n")
  cat(paste(paste("\t", x$run_params[,1], ":", x$run_params[,2]),
            collapse = "\n"))
  cat("\nModel fit statistics:\n")
  cat(paste(paste("\t", x$fit_stats_df[,1], ":", x$fit_stats_df[,2]),
            collapse = "\n"))
  logfile <- !is.null(x$burn_df)
  cat(paste("\nMCMC diagnostics available:", logfile))

}

#' @export
print.structList <- function(x, ...) {
  n <- length(x)
  cat(paste("structList object containing", n, " STRUCTURE runs.\n"))
  cat(paste("Number of Ks by number of runs:"))
  Ks <- table(unlist(lapply(x, getK)))
  print(Ks)
}

#' @export
`[[.structList` <- function(x, ...) {

  y <- NextMethod("[[")
  y
}

#' @export
`[.structList` <- function(x, ...) {
  y <- NextMethod("[")
  # return struct object if list has length 1
  if (length(y) == 1L) {
    return(y)
  } else {
    class(y) <- oldClass(x)
    y
  }
}
