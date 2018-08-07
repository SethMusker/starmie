#' Simulate random Q-matrix and permutations with replications
#'
#' @param n integer number of samples
#' @param alpha vector of length K (number of columns)
#' @param r number of replicates
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats rnorm
#' @export
simulateQ <- function(n, alpha, r) {

  Q_list <- list()
  Q_list[[1]] <- list(Q = rdirichlet(n, alpha),
                      permutations = 1:length(alpha))

  for (i in 2:r) {
    # add noise and shuffle
    q_random <- t(apply(Q_list[[1]]$Q, 1, addNoise))
    size_factor <- rowSums(q_random)
    q_random <- q_random / size_factor
    stopifnot(all.equal(rowSums(q_random), rep(1, n)))
    permutations <- sample(Q_list[[1]]$permutations)

    Q_list[[i]] <- list(Q = q_random[, permutations],
                        permutations = permutations)
  }

  Q_list

}

addNoise <- function(y) {
  1 / (1 + exp(-rnorm(length(y), mean = logit(y),sd = 0.01)))
}

logit <- function(x) log(x / (1-x))
