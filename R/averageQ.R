# averageQ.R

#' Convenience function for creating final Q-matrix after CLUMMP
#' @param Q_list A list of of Q matrices.
#' @export
#' @examples
#' # Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' Q_list <- lapply(multiple_runs_k10, getQ)
#' avgQ <- averageQ(Q_list)
averageQ <- function(Q_list){
  #i/o checks
  if(!all(unlist(lapply(Q_list, inherits, "matrix"))))
    stop("cluster runs must be a list of Q matrices")
  Reduce("+", Q_list) / length(Q_list)
}

