# plotBar.R
# Function for plotting a single Q matrix from Structure or Admixture runs.

#' Generate a barplot of a Structure or Admixture run.
#' @param  x an object of type \code{\link{struct}} or \code{\link{admix}} or a Q-matrix
#' @param  populations a data.frame that contains the sample number as the first column and the population as the second.
#' @param  plot if FALSE returns a data.frame for customised plots
#' @param  facet whether or not to split the barplot by cluster. This is recommended.
#' @import ggplot2
#' @importFrom data.table melt
#' @export
#' @examples
#' # Read file using K = 6 and plot results
#' k6_data <- exampleStructure("barplot")
#' # Generate standard 'structure' barplot
#' plotBar(k6_data, facet = FALSE)
#' # adding group information
#' set.seed(212)
#' pops <- data.frame(id = k6_data$ancest_df[,1],
#' population = sample(letters[1:3], nrow(k6_data$ancest_df), replace = TRUE))
#' # our facetted structure plot
#' plotBar(k6_data, pops)
#' # standard 'structure' bar plot
#' plotBar(k6_data, pops, facet = FALSE)
#' #' admixture example
#' k3_data <- exampleAdmixture()[[3]]
#' plotBar(k3_data)
plotBar <- function(x, populations = NULL, plot = TRUE, facet = TRUE) {
  UseMethod("plotBar", x)
}

#' @method plotBar matrix
#' @export
plotBar.matrix <- function(x, populations = NULL, plot = TRUE, facet = TRUE) {
  #i/o checks
  if (!is.logical(plot)) stop("plot must be one of TRUE or FALSE")
  if (!is.logical(facet)) stop("facet must be one of TRUE or FALSE")

  plot_results <- plotQ(x, populations, facet)
  if (plot) {
    return(plot_results$plot)
  } else {
    return(plot_results$Q_melt)
  }
}

#' @method plotBar struct
#' @export
plotBar.struct <- function(x, populations = NULL, plot = TRUE, facet = TRUE) {
  # for struct if populations = NULL try
  if (is.null(populations)) {
    if ("Pop" %in% colnames(x$ancest_df)) {
      message("Extracting population labels from STRUCTURE output.")
      populations <- x$ancest_df[,c(1,3)]
    }
  }
  Q <- getQ(x)
  plotBar.matrix(Q, populations, plot, facet)
}

#' @method plotBar admix
#' @export
plotBar.admix <- function(x, populations = NULL, plot = TRUE, facet = TRUE ) {
  Q <- getQ(x)
  plotBar.matrix(Q, populations, plot, facet)
}

plotQ <- function(Q, populations_df, facet) {

  if (is.null(populations_df)){
    #Generate a plot without any family information
    Q_melt <- melt(Q, variable.name="Cluster")
    colnames(Q_melt) <- c("Label", "Cluster", "value")
  } else{
    if (!setequal(populations_df[,1], rownames(Q)))
      stop("Mismatch between populations and cluster_run populations.")
    colnames(populations_df) <- c("Label", "Population")
    Q_merge <- merge(populations_df, Q, by.x="Label", by.y=0)
    Q_melt <- melt(Q_merge, id.vars=c("Label", "Population"), variable.name="Cluster")
  }

  #Generate plot
  Q_melt <- Q_melt[order(Q_melt$Cluster),]
  Q_melt$Label <- factor(Q_melt$Label)

  gg <- ggplot(Q_melt, aes_(x=~Label, y=~value, fill=~Cluster))
  if (!is.null(populations_df)){
    if (facet){
      gg <- gg + facet_grid( Cluster ~ Population, scales = "free_x", space = "free_x")
    } else{
      gg <- gg + facet_grid( . ~ Population, scales = "free_x", space = "free_x")
    }
  } else{
    if (facet){
      gg <- gg + facet_grid( Cluster ~ ., scales = "free_x", space = "free_x")
    }
  }
  gg <- gg + geom_bar(stat = "identity", width=1)
  gg <- gg + scale_y_continuous(expand=c(0,0), breaks=c(0.25,0.75))
  gg <- gg + coord_cartesian(ylim=c(0,1))
  gg <- gg + xlab("Sample ID") + ylab("Proportion of cluster")
  gg <- gg + theme_bw()
  gg <- gg + guides(fill=guide_legend(title="Cluster"))
  gg <- gg + theme(axis.text.x = element_text(angle = 90))

  #print plot or return
  return(list(Q_melt = Q_melt, plot=gg))
}
