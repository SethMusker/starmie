# bestK.R
# Functions for determining a suitable K value from multiple Structure runs.
#' Determine a suitable K value from multiple Structure runs
#' @param  x a \code{\link{structList}} or \code{\link{admixList}} object.
#' @param  method the method used to calculate the best K either 'evanno' or
#' 'structure', not required for \code{\link{admixList}} objects.
#' @param  plot whether of not to generate diagnostic plots
#' @import ggplot2
#' @importFrom data.table rbindlist
#' @export
#' @return a data.frame containing  with columns containing the L_k,
#' AIC, BIC, DIC and deltaK for  \code{\link{structList}}. If
#' an \code{\link{admixList}} was given a data.frame returning the log
#' information will be supplied. If plot = TRUE a ggplot object is printed
#' for the method of interest.
#' @details If the K values are not ordered or there an even numbers of runs
#' per K the 'structure' method will be implemented and  the 'evanno' method
#' to compute delta K will not be returned in the output.
#' @examples
#' multi_K <- exampleStructure("multiple_runs")
#' # Run the evanno method and display diagnostic plots.
#' evanno_results <- bestK(multi_K, method = "evanno")
#' # Run the default structure method and display diagnostic plots
#' structure_results <- bestK(multi_K, "structure")
#' # find 'best' K according to results
#' deltaK <- evanno_results$variable == 'delta K'
#' max_deltaK <- which(evanno_results$value == max(evanno_results$value[deltaK], na.rm = TRUE))
#' evanno_results[max_deltaK, ]
#' lK <- structure_results$variable == 'L(K)'
#' max_Lk <- which(structure_results$value == max(structure_results$value[lK], na.rm = TRUE))
#' structure_results[max_Lk,]
#' # admixture example
#' multi_K_admix <- exampleAdmixture()
#' bestK(multi_K_admix)
bestK <- function(x, method, plot = TRUE) {
  UseMethod("bestK", x)
}

#' @method bestK structList
#' @export
bestK.structList <- function(x, method = "evanno", plot = TRUE) {
  if (!(method %in% c("evanno", "structure")))
    stop("method must be one of 'evanno' or 'structure'")
  if ( !is.logical(plot) | is.na(plot) )
    stop("plot must be one of TRUE or FALSE")

  message("Creating diagnostic plots for structure runs.")
  params_sizes <- lapply(x, getD)

  posterior_probs <- data.frame(K=unlist(lapply(x, getK)),
                                pos_prob=unlist(lapply(x, getPosterior)),
                                d = unlist(lapply(params_sizes, function(i) i$d)),
                                n = unlist(lapply(params_sizes, function(i) i$n)))

  # fit stats
  ll_summary <- matrix(unlist(lapply(x, getFitStats)), ncol = 2, byrow = TRUE)
  # compute average AIC/BIC over runs + standard errors
  posterior_probs$aic <- -2*posterior_probs$pos_prob + 2*posterior_probs$d
  posterior_probs$bic <- -2*posterior_probs$pos_prob + posterior_probs$d * log(posterior_probs$n)
  # Gelman 2013 method for computing DIC
  posterior_probs$dic <- -2*ll_summary[,1] + 2*ll_summary[,2]

  # If K are not sequential or do not have equal number of runs just
  # perform the structure approach and produce a simpler plot
  Ks <- table(posterior_probs$K)
  # check for whether it is valid to compute Evanno method
  evanno_invalid <- !is.sequential(names(Ks)) | !all(Ks[[1]]==Ks) | all(Ks == 1)
  if (evanno_invalid | method == "structure") {
    if (method!="structure") {
      if (all(Ks == 1)) {
        stop("Not enough information to compute Evanno statistics.")
      }
      warning("WARNING! K values are not sequential or there are
                an uneven number of runs per K.
                Reverting to structure method instead.")
    }
    model_ll <- bestK_structure(posterior_probs, plot)
    model_ll

  } else {
    model_ll <- bestK_evanno(posterior_probs, plot)
    model_ll
  }
}

#' @method bestK admixList
#' @export
bestK.admixList <- function(x, method = NULL, plot = TRUE) {
  message("Creating diagnositc plots for admixture runs")

  log_df <- combineLogs(x)

  if ( is.null(log_df) ) {
    stop("Need log file information to produce diagnositc plots")
  }
  log_df_tidy <- data.table::melt(combineLogs(x),
                                  id.vars = "K",
                                  variable.name = "statistic")
  if (plot) {
    gg <- ggplot(log_df_tidy, aes_(x = ~K, y = ~value)) +
      geom_point() +
      facet_wrap(~ statistic, ncol = 2, scales = "free_y") +
      theme_bw()
    suppressWarnings(print(gg))
  }
  return(log_df)
}


#' @importFrom stats sd
bestK_evanno <- function(posterior_probs, plot) {
  # plot delta K
  # split by K and calculate the first and second derivatives

  split_posterior_probs <- split(posterior_probs, posterior_probs$K)
  split_posterior_probs <- split_posterior_probs[order(as.numeric(names(split_posterior_probs)))]
  # empty data frame for storing results
  output_variables <- c("L(K)", "L'(K)",  "L''(K)", "delta K")
  output_values <- list()
  k_all <- min(posterior_probs$K):max(posterior_probs$K)
  for (k in k_all) {
    LK <-  split_posterior_probs[[k]]$pos_prob

    if(k > min(k_all) & k < max(k_all)) {
      dK <- split_posterior_probs[[k]]$pos_prob - split_posterior_probs[[k-1]]$pos_prob
      ddK <- split_posterior_probs[[k+1]]$pos_prob-2*split_posterior_probs[[k]]$pos_prob+split_posterior_probs[[k-1]]$pos_prob
      if ( isTRUE(sd(LK) == 0 ) )stop("No deviation between runs of the same K. Evanno statistics cannot be computed.")

    } else {
      dK <- rep(NA, length(LK))
      ddK <- rep(NA, length(LK))
    }


    output_values[[k]] <-  data.frame(K = as.numeric(names(split_posterior_probs)[k]),
                                                variable = output_variables,
                                                value = c(mean(LK), mean(dK), mean(abs(ddK)), mean(abs(ddK))/sd(LK)),
                                                sd = c(sd(LK), sd(dK), sd(ddK), NA))


  }
  posterior_probs_summary <- rbindlist(output_values)

  if (plot) {
    gg <- ggplot(posterior_probs_summary, aes_(x=~K, y=~value)) +
      geom_point() +
      facet_wrap(~variable, ncol=2, scales = "free_y") +
      geom_errorbar(aes_q(ymax=quote(value+sd), ymin=quote(value-sd))) +
      theme_bw() +
      scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K))
    suppressWarnings(print(gg))
  }

  return(posterior_probs_summary)

}

#' @importFrom stats sd
bestK_structure <- function(posterior_probs, plot) {
  # look for change point in log likelihood plot
  # summarise data by K
  posterior_probs_byK <- split(posterior_probs, posterior_probs$K)
  posterior_probs_summary <- rbindlist(lapply(posterior_probs_byK,
                                                     function(i) data.frame(K = unique(i$K),
                                                                            variable = c('L(K)', 'AIC', 'BIC', 'DIC'),
                                                                            value = c(mean(i$pos_prob), mean(i$aic), mean(i$bic), mean(i$dic)),
                                                                            sd = c(sd(i$pos_prob), sd(i$aic), sd(i$bic), sd(i$dic)))))
  #generate plot
  if (plot) {
    gg <- ggplot(posterior_probs_summary, aes_(x=~K, y=~value)) +
      geom_point() +
      geom_errorbar(aes_q(ymax=quote(value+sd),
                        ymin=quote(value-sd)), width=0.1) +
      theme_bw() +
      facet_wrap(~variable, ncol=2, scales = "free_y") +
      scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K)) +
      ylab("Mean log posterior probability")
    suppressWarnings(print(gg))
  }

  return(posterior_probs_summary)

}

is.sequential <- function(x){
  x <- as.numeric(x)
  all(diff(x) == diff(x)[1])
}


