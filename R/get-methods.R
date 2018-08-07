# get-methods.R
# common get-methods for struct and admix objects

#' Retrieve the assumed number of populations from \code{\link{struct}} or
#' \code{\link{admix}} objects.
#' @param x a \code{\link{struct}} object  or \code{\link{admix}} object.
#' @export
getK <- function(x) UseMethod("getK", x)

#' @method getK struct
#' @export
getK.struct <- function(x) {
  x$K
}

#' @method getK admix
#' @export
getK.admix <- function(x) {
  x$K
}

#' @method getK fast
#' @export
getK.fast <- function(x) {
  x$K
}

#' Retrieve Q matrix from \code{\link{struct}} or \code{\link{admix}} \code{\link{fast}} objects.
#'
#' @param x a \code{\link{struct}} or \code{\link{admix}} or \code{\link{fast}} object.
#' @export
getQ <- function(x) {
  UseMethod("getQ", x)
}

#' @method getQ struct
#' @export
getQ.struct <- function(x) {
  columns_to_keep <- colnames(x$ancest_df)
  columns_to_keep <- grep("^Cluster", columns_to_keep)
  Q <- data.matrix(x$ancest_df[, columns_to_keep])
  rownames(Q) <- x$ancest_df$Label
  return(Q)
}

#' @method getQ admix
#' @export
getQ.admix <- function(x) {
  as.matrix(x$Q_df)
}

#' @method getQ fast
#' @export
getQ.fast <- function(x) {
  as.matrix(x$Q_df)
}

#' Retrieve estimated within-cluster allele frequencies
#'
#' @param x a \code{\link{struct}} or \code{\link{admix}} or \code{\link{fast}} object.
#' @export
getClusterAlleleFreqMat <- function(x) {
  UseMethod("getClusterAlleleFreqMat", x)
}

#' @method getClusterAlleleFreqMat struct
#' @export
getClusterAlleleFreqMat.struct <- function(x) {
  do.call("rbind",
          lapply(x$clust_allele_list,
                 function(x) cbind(rep(x$Locus, x$AlleleNumber), x$FreqMatrix[, -2])))

}

#' @method getClusterAlleleFreqMat admix
#' @export
getClusterAlleleFreqMat.admix <- function(x) {
  x$P_df
}

#' @method getClusterAlleleFreqMat fast
#' @export
getClusterAlleleFreqMat.fast <- function(x) {
  x$P_df
}

#' Retrieve estimated population allele frequencies
#'
#' @param x a \code{\link{struct}} object.
#' @importFrom data.table rbindlist
#' @export
getCompleteAlleleFreqMat <- function(x) {
  UseMethod("getCompleteAlleleFreqMat", x)
}
getCompleteAlleleFreqMat.struct <- function(x) {
  do.call("rbind", lapply(x$clust_allele_list,
                                function(x) cbind(rep(x$Locus, x$AlleleNumber),
                                                  x$FreqMatrix[, c(1,2)])))
}
