# plot the pairwise coancestry matrix

#' fineSTRUCTURE pairwise coancestry plot
#'
#' @param x a \code{\link{fineStruct}} object
#' @param cex (default=0.7) number indicating the amount by which plotting text and symbols should be scaled relative to the default.
#' @param colours a vector of colors to be passed to image. (e.g., heat.colors(n=200)[200:1])
#' @param keep_ones (default=FALSE) whether or not to annotate node with full support.
#' @param max_coancestry the maximum value to be reprsented in the heatmap. Large values will be reduced to max_coancestry.
#' @export
#' @examples
#' # fineStruct example
#' chunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.chunkcounts.out", package = "starmie")
#' treefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.tree.xml", package = "starmie")
#' mcmcfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mcmc.xml", package = "starmie")
#' meancoincidencefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.meancoincidence.csv", package = "starmie")
#' mappopchunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mapstate.csv", package = "starmie")
#' fineData <- loadFineStructure(chunkfile, treefile, mcmcfile, meancoincidencefile, mappopchunkfile)
#' plotPWCoancestry(fineData, max_coancestry=500)
plotPWCoancestry <- function(x, cex=0.5, colours=NULL, keep_ones=FALSE, max_coancestry=NULL) {
  #i/o checks
  if(!(inherits(x, "fineStruct"))) stop("x must be a fineStruct object")
  if(!is.numeric(cex)) stop("cex must be numeric")
  if(!is.logical(keep_ones)) stop("keep_ones must be one of TRUE or FALSE")
  if(!is.numeric(max_coancestry) & !is.null(max_coancestry)) stop("max_coancestry must be numeric")
  if(is.null(x$chunkcounts_matrix)) stop("fineStruct object has no chunkcounts_matrix")
  if(!requireNamespace("phytools", quietly = TRUE)) stop("phytools package not installed, please install it")

  if (is.null(colours)){
    colours <- grDevices::colorRampPalette(colors=c("#e0ecf4","#8856a7"))(200)
  }

  #match orderings
  t_matrix <- x$chunkcounts_matrix

  if (!is.null(max_coancestry)){
    t_matrix[t_matrix>max_coancestry] <- max_coancestry
  }

  diag(t_matrix) <- NA
  t_matrix <- t_matrix[match(x$dendro$tip.label,rownames(t_matrix)),]
  t_matrix <- t_matrix[,match(x$dendro$tip.label,colnames(t_matrix))]

  node_labels <- x$dendro$node.label
  if (!keep_ones){
    node_labels[as.numeric(node_labels)==1] <- ""
  }

  phytools::phylo.heatmap(x$dendro, t_matrix, fsize=cex
                          , colors=colours)
  ape::nodelabels(node_labels, node=1:x$dendro$Nnode+ape::Ntip(x$dendro),
                  adj=c(1,-0.2),frame="none", cex=cex)

}




