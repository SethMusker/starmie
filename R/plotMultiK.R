# plotMultiK.R
# Function for plotting Q matrices for multiple K values from Structure or Admixture runs.

#' Generate a barplot for multiple values of K..
#' @param  x a \code{\link{structList}} or \code{\link{admixList}} object or a list of Q-matrices.
#' @param  populations a data.frame that contains the sample number as the first column and the population as the second.
#' @param  plot if FALSE returns a data.frame for customised plots
#' @import ggplot2
#' @importFrom data.table melt
#' @export
#' @examples
#' cluster_runs <- exampleStructure("multiple_runs")
#' # Generate barplot
#' plotMultiK(cluster_runs[3:5])
plotMultiK <- function(x, populations = NULL, plot = TRUE) {
  UseMethod("plotMultiK", x)
}

#' @method plotMultiK list
#' @export
plotMultiK.list <- function(x, populations = NULL, plot = TRUE) {
  # i/o checks
  if (!all(unlist(lapply(x, inherits, "matrix"))) ) {
    stop("All elements of list must be matrix objects.")
  }
  plot_results <- plotMultiQ(x, populations)
  if (plot) {
    return(plot_results$plot)
  } else {
    return(plot_results$Q_melt)
  }
}

#' @method plotMultiK structList
#' @export
plotMultiK.structList <- function(x, populations = NULL, plot = TRUE ) {
  Q_list <- lapply(x, getQ)

  #If no population labels are given try to use
  # those available in the structure run.
  if (is.null(populations)) {
    if ("Pop" %in% colnames(x[[1]]$ancest_df)) {
      message("Extracting population labels from STRUCTURE output.")
      populations <- x[[1]]$ancest_df[,c(1,3)]
    }
  }
  plotMultiK.list(Q_list, populations, plot)

}


#' @method plotMultiK admixList
#' @export
plotMultiK.admixList <- function(x, populations, plot = TRUE) {
  Q_list <- lapply(x, getQ)
  plotMultiK.list(Q_list, populations, plot)

}

#' @importFrom stats ave
plotMultiQ <- function(Q_list, populations_df){

  #get K labels
  Ks <- unlist(lapply(Q_list, ncol))
  if (length(unique(Ks))!= length(Ks)){
    #Repeated Ks so label with subnumbering
    Ks <- paste(Ks, ave(Ks, Ks, FUN=seq_along), sep=".")
  } else{
    Ks <- as.character(Ks)
  }

  for (i in 1:length(Q_list)){
    if (is.null(rownames(Q_list[[i]]))) {
      Label <- 1:nrow(Q_list[[i]])
    } else {
      Label <- rownames(Q_list[[i]])
    }

    Q_list[[i]] <- data.frame(Label=Label, Q_list[[i]],
                              K=rep(Ks[[i]], nrow(Q_list[[i]])))
  }

  #Melt and append Q matrices
  Q_melt <- do.call("rbind", lapply(Q_list, melt, id.vars=c("Label","K"), variable.name="Cluster"))

  if (!is.null(populations_df)){
    #Generate a plot with family information
    colnames(populations_df) <- c("Label", "Population")
    Q_melt <- merge(populations_df, Q_melt, by.x="Label", by.y="Label")
  }

  #Generate plot
  Q_melt <- Q_melt[order(Q_melt$Cluster),]
  Q_melt$Label <- factor(Q_melt$Label)
  gg <- ggplot(Q_melt, aes_(x=~Label, y=~value, fill=~Cluster))
  if (!is.null(populations_df)){
    gg <- gg + facet_grid( K ~ Population, scales = "free_x", space = "free_x")
  } else{
    gg <- gg + facet_grid( K ~ ., scales = "free_x", space = "free_x")
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
