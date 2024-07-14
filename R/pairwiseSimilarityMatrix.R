# pairwiseSimilarityMatrix.R

#' Compute matrix of pairwise similarity between Q matrices.
#' @param Q_list A list of of Q matrices.
#' @details Implementation of the pairwise similarity metric as
#' defined in Jakobsson, M. and Rosenberg, N. A., 2007.
#' @export
#' @examples
#' # Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' Q_list <- lapply(multiple_runs_k10, getQ)
#' H <- pairwiseSimilarityMatrix(Q_list)
pairwiseSimilarityMatrix<-function(Q_list){
  R <- length(Q_list)
  combs <- combn(1:R,2)
  H <- matrix(nrow=R,ncol=R)
  for(i in 1:ncol(combs)){
    a<-combs[1,i]
    b<-combs[2,i]
    H[a,b]<-G(Q_list[[a]],Q_list[[b]])
    H[b,a]<-H[a,b]
  }
  diag(H)<-1
  return(H)
}
