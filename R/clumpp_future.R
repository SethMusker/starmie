clumpp_future<-function (Q_list, method = "greedy", iter = 100) 
{
  require(future.apply)
  plan(multisession)
  require(tictoc)
  tic()
  if (!(method %in% c("greedy", "greedyLargeK", "stephens"))) {
    stop("Not a valid CLUMPP method, please use on of: 'greedy', 'greedyLargeK' or 'stephens'")
  }
  if (!all(unlist(lapply(Q_list, inherits, "matrix")))) 
    stop("cluster runs must be a list of Q matrices")
  if (!all.equal(iter, as.integer(iter)) || iter < 0) 
    stop("number of iterations must be a positive integer")
  if (length(Q_list) == 1) {
    message("Q_list only contains 1 element, clumpping is not necessary.")
    return(Q_list)
  }
  dim_Q <- dim(Q_list[[1]])
  if (!all(unlist(lapply(Q_list[-1], function(x) all(dim(x) == 
                                                     dim_Q))) == TRUE)) 
    stop("size of all matrices in Q_list must be equal")
  K <- dim_Q[2]
  if (K == 1) {
    message("Number of assumed populations = 1 for all Q matrices, clummping is unecessary, returning Q_list")
    return(Q_list)
  }
  if (method == "greedy") {
    perms <- replicate(iter, sample(1:length(Q_list), size = length(Q_list), 
                                    replace = FALSE))
    if (K > 8) {
      permQs <- future_apply(perms, 2, function(p) iterativeGreedy(Q_list[p]))
      Hs <- future_lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    }
    else {
      permQs <- future_apply(perms, 2, function(p) memoryGreedy(Q_list[p]))
      Hs <- future_lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    }
  }
  else if (method == "greedyLargeK") {
    perms <- replicate(iter, sample(1:length(Q_list), size = length(Q_list), 
                                    replace = FALSE))
    permQs <- future_apply(perms, 2, function(p) largeKGreedy(Q_list[p]))
    Hs <- future_lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
    Q_list <- permQs[[which.max(Hs)]]
  }
  else if (method == "stephens") {
    Q_list <- getStephens(Q_list)
  }
  return(Q_list)
  toc()
}
environment(my_clumpp) <- asNamespace('starmie')
