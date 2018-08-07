# plotTreeBar.R
# Function for plotting a single Q matrix with a dendrogram from Structure or Admixture runs.

#' Generate a barplot of a Structure or Admixture run.
#' @param  x a single cluster run object of type \code{\link{struct}} or \code{\link{admix}} or Q-matrix
#' @param  facet whether or not to split the barplot by cluster. This is recommended.
#' @param  dendro an object of class 'hclust' (defaults to hclust with average linkage)
#' @param  cut an integer vector output by 'cutree' (defaults to cutree k=ncols(Q))
#' @import ggplot2
#' @import gridExtra
#' @import ggdendro
#' @importFrom data.table melt
#' @export
#' @examples
#' # Read file using K = 6 and plot results
#' k6_data <- exampleStructure("barplot")
#' # our facetted structure plot with tree
#' plotTreeBar(k6_data)
#' # standard 'structure' bar plot with tree
#' plotTreeBar(k6_data, facet = FALSE)
#' # Admix example
#' k3_data <- exampleAdmixture()[[3]]
#' plotTreeBar(k3_data)
plotTreeBar <- function(x, facet = TRUE, dendro = NULL, cut = NULL ) {
  UseMethod("plotTreeBar", x)
}

#' @method plotTreeBar matrix
#' @importFrom stats hclust cutree
#' @export
plotTreeBar.matrix <- function(x, facet = TRUE, dendro = NULL, cut = NULL) {
  #If no clustering is given default to average linkage
  if (!inherits(dendro, "hclust") & !is.null(dendro))
    stop("dendro must be a hclust object")

  if (!inherits(cut, "integer") & !is.null(cut))
    stop("cut must be an integer vector")

  if (is.null(dendro)) {
    dendro <- hclust(dist(x, method="euclidean")
                     , method="average")
  }
  if (is.null(cut)) {
    cut <- cutree(dendro, k=ncol(x))
  }

  plot_results <- plotTreeQ(x, facet, dendro, cut)

  return(plot_results)
}

#' @method plotTreeBar struct
#' @export
plotTreeBar.struct <- function(x, facet = TRUE, dendro = NULL, cut = NULL) {
  Q <- getQ(x)
  plotTreeBar.matrix(Q, facet, dendro, cut)
}

#' @method plotTreeBar admix
#' @export
plotTreeBar.admix <- function(x, facet = TRUE, dendro = NULL, cut = NULL) {
  Q <- getQ(x)
  plotTreeBar.matrix(Q, facet, dendro, cut)
}

plotTreeQ <- function(Q, facet, dendro, cut) {

  # if dendro labels is null set to row ordering of Q
  if (is.null(rownames(Q))) rownames(Q) <- as.character(1:nrow(Q))
  if( is.null(dendro$labels) ) dendro$labels <- as.character(1:nrow(Q))
  populations_df <- data.frame(Label=factor(dendro$labels,
                                            levels=dendro$labels[dendro$order]),
                               Population=cut)

  if (!setequal(populations_df[,1], rownames(Q)))
    stop("Mismatch between Q and dendrogram labels")

  colnames(populations_df) <- c("Label", "Population")

  #Generate plot

  ##Dendrogram
  row.dendro <- dendro_data(dendro, type="rectangle")
  row.plot <- mydplot(row.dendro, col=TRUE, labels=TRUE) +
    scale_x_continuous(breaks = 1:length(dendro$labels),
                       labels=dendro$labels[dendro$order]) +
    theme(plot.margin = unit(c(0,0,0,0), "lines"),
          axis.ticks = element_blank(), axis.text.x = element_blank())
  row.ord <- match(dendro$labels[dendro$order], rownames(Q))

  Q_merge <- merge(populations_df, Q, by.x="Label", by.y=0)
  Q_merge <- merge(Q_merge, row.dendro$labels[,-2], by.x="Label", by.y='label')
  Q_melt <- melt(Q_merge, id.vars=c("Label", "Population", "x"),
                 variable.name="Cluster")

  ##Barplot
  gg <- ggplot(Q_melt, aes_(x=~x, y=~value, fill=~Cluster))

  if (facet){
    gg <- gg + facet_grid( Cluster ~ ., scales = "free_x", space = "free_x")
  }

  gg <- gg + geom_bar(stat = "identity", width=1)
  gg <- gg + scale_y_continuous(expand=c(0,0), breaks=c(0.25,0.75))
  gg <- gg + scale_x_continuous(breaks = 1:length(dendro$labels),
                                labels=dendro$labels[dendro$order])
  gg <- gg + coord_cartesian(ylim=c(0,1))
  gg <- gg + xlab("Sample ID") + ylab("Proportion of cluster")
  gg <- gg + theme_bw()
  gg <- gg + guides(fill=guide_legend(title="Cluster"))
  gg <- gg + theme(axis.text.x = element_text(angle = 90),
                   plot.margin = unit(c(0,0,0,0), "lines")
                   )
  change_points <- populations_df[match(dendro$labels[dendro$order],
                                        populations_df$Label),]$Population
  change_points <- c(c(1,1+which(diff(change_points)!=0)),
                     length(dendro$labels)+1)-0.5

  gg <- gg + geom_vline(xintercept = change_points)
  if (facet){
    gg <- gg + theme(legend.position="none")
  }

  ##Merge plots
  gp1 <- ggplotGrob(row.plot)
  gp2 <- ggplotGrob(gg)
  maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
  gp1$widths[2:5] <- as.list(maxWidth)
  gp2$widths[2:5] <- as.list(maxWidth)

  gg <- grid.arrange(gp1, gp2, ncol=1,heights=c(2/5,3/5))

  #print plot or return
  return(gg)
}


##Taken from Chris Wallace's heatmap code
##https://gist.github.com/chr1swallace/4672065
mydplot <- function(ddata, row=!col, col=!row, labels=col) {
  ## plot a dendrogram
  yrange <- range(ddata$segments$y)
  yd <- yrange[2] - yrange[1]
  nc <- max(nchar(as.character(ddata$labels$label)))
  tangle <- if(row) { 0 } else { 90 }
  tshow <- col
  p <- ggplot() +
    geom_segment(data=segment(ddata), aes_(x=~x, y=~y, xend=~xend, yend=~yend)) +
    labs(x = NULL, y = NULL) + theme_dendro()
  if(row) {
    p <- p +
      scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
      coord_flip()
  } else {
    p <- p +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  return(p)
}
