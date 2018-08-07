# plot-methods for visualising full Structure model

#' Plot principal coordinates from Q-matrix, struct or admix objects
#'
#' @param x a Q-matrix of probability memberships, or \code{\link{struct}} or \code{\link{admix}} object
#' @param method (default = NULL) string either 'nnd' or 'jsd' valid only for \code{\link{struct}} objects
#' @details "nnd" uses the nucleotide distance matrix estimated by STRUCTURE
#' to construct the principal coordinates, sizing the points by the expected
#' heterozygosity within a cluster. "jsd" produces a principal coordinates
#' from the Jensen Shannon Divergence metric as used by the 'ldavis' package and
#' is the default for Q-matrix or admix objects. By default using plotMDS on
#' a struct object will produce principal coordinates on the clusters
#' themselves rather than within samples.
#' @importFrom proxy dist
#' @import ggplot2
#' @importFrom stats cmdscale
#' @export
#' @examples
#' # struct example
#' k6_data <- exampleStructure("barplot")
#' plotMDS(k6_data)
#' plotMDS(k6_data, method = "jsd")
#' # admix example
#' k3_data <- exampleAdmixture()[[3]]
#' plotMDS(k3_data)
plotMDS <- function(x, method = NULL) {
    UseMethod("plotMDS", x)
}

#' @method plotMDS matrix
#' @export
plotMDS.matrix <- function(x, method = NULL) {

  dist_xy <- proxy::dist(x, method = .JSD)
  mds_clust <- cmdscale(dist_xy)
  mds_df <- data.frame(PC1 = mds_clust[,1],
                       PC2 = mds_clust[,2])
  ggplot(mds_df, aes_(x = ~PC1, y = ~PC2)) +
    geom_point() +
    theme_bw()
}

#' @method plotMDS struct
#' @export
plotMDS.struct <- function(x, method = "nnd") {
  # gather visual elements
  if (!(method %in% c("nnd", "jsd")))
    stop("Not a valid method, must be either 'nnd' or 'jsd'")
  if (method == "nnd") {
    # use structure allele-freqs diveragnes as distance matrix
    dist_xy <- x$allele_freqs
    # set diags to 0
    diag(dist_xy) <- 0
    # grab cluster expected heterzygosities
    clust_eh <- x$avg_dist_df

    # compute mds coordinates
    mds_clust <- cmdscale(dist_xy)

    mds_df <- data.frame(clust_eh, PC1 = mds_clust[,1], PC2 = mds_clust[,2])

    plot_out <- ggplot(mds_df, aes_(x = ~PC1, y = ~PC2, size = ~Avg.dist)) +
      geom_point() +
      guides(size = guide_legend(title = "Expected Heterozygosity"))+
      theme_bw()
  }
  if (method == "jsd") {
    Q <- getQ(x)
    dist_xy <- proxy::dist(t(Q), method = .JSD)
    mds_clust <- cmdscale(dist_xy)
    mds_df <- data.frame(Cluster = as.integer(sub("Cluster ", "", colnames(Q))),
                         Relative.Contribution = colSums(Q) / sum(Q),
                         PC1 = mds_clust[,1],
                         PC2 = mds_clust[,2])
    plot_out <- ggplot(mds_df,
                       aes_(x = ~PC1, y = ~PC2, size = ~Relative.Contribution)) +
      geom_point() +
      guides(size = guide_legend(title = "Relative Cluster Contribution")) +
      theme_bw()
  }
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot_out + ggrepel::geom_text_repel(aes_(label = ~Cluster),
                                        size = 4,
                                        box.padding = unit(0.35, "lines"),
                                        point.padding = unit(0.3, "lines"))
  } else {
    plot_out
  }
}

#' @method plotMDS admix
#' @export
plotMDS.admix <- function(x, method = NULL) {
  Q_hat <- getQ(x)
  plotMDS.matrix(Q_hat)
}

.JSD<- function(x,y) sqrt(0.5 * .KLD(x, (x+y)/2) + 0.5 * .KLD(y, (x+y)/2))
.KLD <- function(x,y) sum(x * log(x/y))
