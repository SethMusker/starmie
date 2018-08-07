# plot the Pairwise coincidence, showing the MAP state

#' fineSTRUCTURE pairwise coincidence plot
#'
#' @param x a \code{\link{fineStruct}} object
#' @param cex (default=0.7) number indicating the amount by which plotting text and symbols should be scaled relative to the default.
#' @param keep_ones (default=FALSE) whether or not to annotate node with full support.
#' @param colours a vector of colors to be passed to image. (e.g., heat.colors(n=200)[200:1])
#' @export
#' @examples
#' # fineStruct example
#' chunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.chunkcounts.out", package = "starmie")
#' treefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.tree.xml", package = "starmie")
#' mcmcfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mcmc.xml", package = "starmie")
#' meancoincidencefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.meancoincidence.csv", package = "starmie")
#' mappopchunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mapstate.csv", package = "starmie")
#' fineData <- loadFineStructure(chunkfile, treefile, mcmcfile, meancoincidencefile, mappopchunkfile)
#' plotPWCoincidence(fineData)
plotPWCoincidence <- function(x, cex=0.5, keep_ones=FALSE, colours=NULL) {
  #i/o checks
  if(!(inherits(x, "fineStruct"))) stop("x must be a fineStruct object")
  if(!is.numeric(cex)) stop("cex must be numeric")
  if(!is.logical(keep_ones)) stop("keep_ones must be one of TRUE or FALSE")
  if(is.null(x$meancoincedence_matrix)) stop("fineStruct object has no meancoincedence_matrix")
  if(!requireNamespace("phytools", quietly = TRUE)) stop("phytools package not installed, please install it")

  if (is.null(colours)){
    colours <- grDevices::colorRampPalette(colors=c("blue","red"))(100)
  }

  #match orderings
  t_matrix <- x$meancoincedence_matrix
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




