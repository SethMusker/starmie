# plots a basic dendrogram of the fineStructure ouput

#' Plot dendrogram from output of fineStructure
#'
#' @param x a \code{\link{fineStruct}} or ape::phylo object
#' @param cex (default=0.7) number indicating the amount by which plotting text and symbols should be scaled relative to the default.
#' @param ... Additional arguments to be passed to ape::plot.phylo
#' @export
#' @examples
#' # fineStruct example
#' chunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.chunkcounts.out", package = "starmie")
#' treefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.tree.xml", package = "starmie")
#' mcmcfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mcmc.xml", package = "starmie")
#' fineData <- loadFineStructure(chunkfile, treefile, mcmcfile)
#' plotFineStructDendro(fineData)
plotFineStructDendro <- function(x, cex=0.7, ...) {
  #i/o checks
  if(!requireNamespace("ape", quietly = TRUE)) stop("ape package not installed, please install it")

  UseMethod("plotFineStructDendro", x)
}

#' @method plotFineStructDendro phylo
#' @export
plotFineStructDendro.phylo <- function(x, cex=0.7, ...) {
  ape::plot.phylo(x, direction = "right", label.offset = 0.1, cex = cex, ...)
  ape::nodelabels(x$node.label, frame = "none", adj = c(2,-0.5), cex=cex)
}

#' @method plotFineStructDendro fineStruct
#' @export
plotFineStructDendro.fineStruct <- function(x, cex=0.7, ...) {
  plotFineStructDendro.phylo(x$dendro, cex=cex, ...)
}




