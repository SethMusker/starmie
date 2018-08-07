# plotMCMC.R
# Function for plotting MCMC chains for diagnostics

#' Plot STRUCTURE MCMC chains
#'
#' @param x \code{\link{structList}} objects or a \code{\link{struct}} object
#' @param plot logical print resulting plot default TRUE
#' @param facet logical facet by K default TRUE
#' @param use_logL logical plot log-likelihood (TRUE) or admixture coefficient
#' @description Plot non-burn MCMC iterations of STRUCTURE for checking convergence.
#' If plot is set to FALSE a data.frame is returned containing the log likelihood
#' and alpha values over different K and runs and not plot is printed to the device.
#' @return If plot is TRUE a ggplot is printed to the screen and the
#' plot object and the data to generate it are returned. Otherwise,
#' a data.frame containing MCMC info it returned.
#' @import ggplot2
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' #Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' print(multiple_runs_k10)
#' results <- plotMCMC(multiple_runs_k10, plot = TRUE)
#' single_run <- plotMCMC(multiple_runs_k10[[1]])
plotMCMC <- function(x, plot = TRUE, use_logL = TRUE, facet = TRUE) {
  UseMethod("plotMCMC", x)
}

#' @method plotMCMC struct
#' @export
plotMCMC.struct <- function(x, plot = TRUE, use_logL = TRUE, facet = NULL) {
  stopifnot(is.logical(plot))
  stopifnot(is.logical(use_logL))
  mcmc_df <- getMCMC(x)

  if (plot) {
    if (use_logL) {
      gg <- ggplot(mcmc_df,
                   aes_(x = ~Iteration, y = ~LogL))
    } else {
      gg <- ggplot(mcmc_df,
                   aes_(x = ~Iteration, y = ~Alpha))
    }

    gg <-  gg +
      geom_line() + theme_bw()

    print(gg)
    return(list(mcmc_info = mcmc_df, mcmc_plot = gg))
  }
  mcmc_df

}

#' @method plotMCMC structList
#' @importFrom stats ave
#' @export
plotMCMC.structList <- function(x, plot = TRUE, use_logL = TRUE, facet = TRUE) {


  # generate data frame of mcmc diagnostics from non-burn iterations
  mcmc_df <- rbindlist(lapply(x, getMCMC))
  # run number is kind of arbirtary here, just group by K and iteration
  mcmc_df$run <- factor(ave(mcmc_df$K, mcmc_df$Iteration, FUN = seq_along))

  if (plot) {
    if (use_logL) {
      gg <- ggplot(mcmc_df,
                   aes_(x = ~Iteration, y = ~LogL, colour = ~run))
    } else {
      gg <- ggplot(mcmc_df,
                   aes_(x = ~Iteration, y = ~Alpha, colour = ~run))
    }

    gg <-  gg +
      geom_line() +
      scale_colour_hue("Run No.") +
      theme(legend.position = "bottom") + theme_bw()

    if (facet) {
      gg <- gg + facet_grid(K ~ . , scales = "free_y")
    }

    print(gg)
    return(list(mcmc_info = mcmc_df, mcmc_plot = gg))
  }
  mcmc_df
}
