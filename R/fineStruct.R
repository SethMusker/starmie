# fineStruct.R
#' Constructor for fineStruct object
#'
#' @return an fineStruct object which is a list with 3 elements:
#'  nsamples: number of samples used
#'
#'  cfactor: estimate of the c parameter in the fineStructure model
#'
#'  chunkcounts_matrix: a data.matrix of chunkcounts
#'
#'  dendro: a dendrogram object estimated by fineStructure
#'
#'  mcmc_df: a data.frame containing the mcmc output of fineStructure
#'
#'  meancoincedence_matrix: a data.matrix of pairwise coincidence, showing the MAP state
#'
#'  mappopchunk_matrix: a data.matrix of population-by-population chunkcounts
#'
#' @export
fineStruct <- function() {
  structure(list(nsamples = NULL, chunkcounts_matrix = NULL, cfactor = NULL,
                 dendro =  NULL, mcmc_df = NULL, meancoincedence_matrix=NULL,
                 mappopchunk_matrix=NULL),
            class = "fineStruct")
}

#' @export
print.fineStruct <- function(x, ...) {
  cat("fineStruct object containing run information\n")
  cat(paste("  ", "No. samples:", x$nsamples, "\n"))
  cat("Chunk count dataframe:\n")
  print(head(x$chunkcounts_matrix, n=2))
  cat("Dendrogram:\n")
  print(x$dendro)
  cat("MCMC dataframe:\n")
  print(head(x$mcmc_df, n=2))
  if (!is.null(x$meancoincedence_matrix)) print("Includes mean coincedence matrix")
  if (!is.null(x$mappopchunk_matrix)) print("Includes population-by-population chunkcount matrix")
}
