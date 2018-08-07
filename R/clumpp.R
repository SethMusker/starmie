# clumpp.R
# Functions for merging Q matrices from Structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The algorithm to use to infer the correct permutations. One of 'greedy' or 'greedyLargeK' or 'stephens'
#' @param iter The number of iterations to use if running either 'greedy' or 'greedyLargeK'
#' @import iterpc
#' @importFrom combinat permn
#' @importFrom purrr map_dbl
#' @export
#' @examples
#' # use multiple K=3 runs
#' cl_data <- exampleStructure("clumpp")
#' print(cl_data)
#' Q_list <- lapply(cl_data, getQ)
#' clumppy <- clumpp(Q_list)
clumpp <- function(Q_list, method="greedy", iter=100){

  # i/o checks
  if (!(method %in% c("greedy", "greedyLargeK", "stephens"))) {
    stop("Not a valid CLUMPP method, please use on of: 'greedy', 'greedyLargeK' or 'stephens'")
  }

  if(!all(unlist(lapply(Q_list, inherits, "matrix"))))
    stop("cluster runs must be a list of Q matrices")

  if(!all.equal(iter, as.integer(iter)) || iter<0) stop("number of iterations must be a positive integer")

  # if length of Q_list is 1, clummping is not necessary
  if (length(Q_list) == 1) {
    message("Q_list only contains 1 element, clumpping is not necessary.")
    return(Q_list)
  }

  # check dims of input Q_list are all equal
  dim_Q <- dim(Q_list[[1]])
  if(!all(unlist(lapply(Q_list[-1], function(x) all(dim(x) == dim_Q))) == TRUE))
    stop("size of all matrices in Q_list must be equal")

  K <- dim_Q[2]

  # gracefully return Q_list if K=1
  if(K == 1) {
    message("Number of assumed populations = 1 for all Q matrices, clummping is unecessary, returning Q_list")
    return(Q_list)
  }

  if (method=="greedy"){
    #Greedy clumpp algorithm
    perms <- replicate(iter, sample(1:length(Q_list),size=length(Q_list),replace=FALSE))
    if (K>8){
      permQs <- apply(perms, 2, function(p) iterativeGreedy(Q_list[p]))
      Hs <- lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    }
    else{
      permQs <- apply(perms, 2, function(p) memoryGreedy(Q_list[p]))
      Hs <- lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
      Q_list <- permQs[[which.max(Hs)]]
    }

  } else if (method=="greedyLargeK"){
    #Use LargeKGreedy algorithm
    perms <- replicate(iter, sample(1:length(Q_list),size=length(Q_list),replace=FALSE))
    permQs <- apply(perms, 2, function(p) largeKGreedy(Q_list[p]))
    Hs <- lapply(permQs, function(x) averagePairWiseSimilarityH(x$Q_list))
    Q_list <- permQs[[which.max(Hs)]]

  } else if (method == "stephens") {
    Q_list <- getStephens(Q_list)

  }
  return(Q_list)
}

#' Use the Stephen's method to permute sample labels
#' @param Q_list A list of of Q matrices.
getStephens <- function(Q_list){
  # Create 3-dimensional array for input into stephens function
  # dimensions are equal to R by n by K
  # R := number of runs
  # n := number of rows in Q matrix
  # K := number of columns in Q matrix
  # convert list to array then transpose columms
  p <- aperm(simplify2array(Q_list), c(3,1,2))

  if(requireNamespace("label.switching", quietly = TRUE)) {
    perm <- label.switching::stephens(p)
  } else {
    stop("label.switching package not installed, please install it")
  }

  # reorder columns in according to new permuations
  # Rename columns
  column_names <- paste("Cluster ", seq_len(dim(p)[3]))
  Q_update <- lapply(seq_len(dim(p)[1]),
         function(i) {
           q_perm <- Q_list[[i]][, perm$permutations[i, ]]
           colnames(q_perm) <- column_names
           q_perm
         }
  )


  return(list(Q_list=Q_update, permutations=perm$permutations))
}

memoryGreedy <- function(Q_list){
  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow=length(Q_list), ncol=K)
  permutations[1,] <- seq(1,K)

  #Faster but with a high memory footprint for large K
  for (i in 1:(length(Q_list)-1)){
    permuations <- permn(1:ncol(Q_list[[i+1]]))
    perm_scores <- map_dbl(permuations, J_perm, Q_list[[i+1]], Q_list[1:i])
    perm <- permuations[[which.max(perm_scores)]]
    Q_list[[i+1]] <- Q_list[[i+1]][,perm]
    permutations[i+1,] <- perm
  }

  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x){ length(unique(colnames(x)))==ncol(x) })))) stop("Duplicated column names in output Q matrices")

  #Rename columns
  column_names <- paste("Cluster ", seq(1,ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x){
    colnames(x) <- column_names
    return(x)})

  return(list(Q_list=Q_list, permutations=permutations))
}

iterativeGreedy <- function(Q_list){
  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow=length(Q_list), ncol=K)
  permutations[1,] <- seq(1,K)

  #slower but memory efficient
  for (i in 1:(length(Q_list)-1)){
    permuations <- iterpc(ncol(Q_list[[i+1]]), ordered = TRUE)

    max_perm <- 1:ncol(Q_list[[i+1]])
    max <- -Inf
    j <- 0
    while(j<getlength(permuations)) {
      value <- J_perm(getnext(permuations), Q_list[[i+1]], Q_list[1:i])
      if (value>max){
        max <- value
        max_perm <- getcurrent(permuations)
      }
      j <- j+1
    }
    Q_list[[i+1]] <- Q_list[[i+1]][,max_perm]
    permutations[i+1,] <- max_perm
  }

  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x){ length(unique(colnames(x)))==ncol(x) })))) stop("Duplicated column names in output Q matrices")

  #Rename columns
  column_names <- paste("Cluster ", seq(1,ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x){
    colnames(x) <- column_names
    return(x)})

  return(list(Q_list=Q_list, permutations=permutations))
}

largeKGreedy <- function(Q_list){
  #Initial iteration
  #calculate pairwise column comparisons
  column_pairs <- getall(iterpc(ncol(Q_list[[1]]), 2, ordered = TRUE, replace = TRUE))
  pair_comparisons <- apply(column_pairs, 1
                            , function(x) { list(G = G(Q_list[[1]][,x[1],drop=FALSE], Q_list[[2]][,x[2],drop=FALSE])
                                                 , Qy = x[1]
                                                 , Qz = x[2])})
  pair_comparisons_df <- data.frame(matrix(unlist(pair_comparisons), ncol=3, byrow=TRUE))

  #Create matrix to store permutations
  K <- ncol(Q_list[[1]])
  permutations <- matrix(0, nrow=length(Q_list), ncol=K)
  permutations[1,] <- seq(1,K)

  #permute the second Q matrix
  perm <- get_best_permutation(pair_comparisons_df)
  Q_list[[2]] <- Q_list[[2]][,perm]
  permutations[2,] <- perm

  if (length(Q_list)>2){
    for (i in 3:length(Q_list)){
      #remaining iterations
      pair_comparisons <- apply(column_pairs, 1
                                , function(x) { list(J = J_largeK(Q_list[1:(i-1)], Q_list[[i]], x[1], x[2])
                                                     , Qy = x[1]
                                                     , Qz = x[2])})
      pair_comparisons_df <- data.frame(matrix(unlist(pair_comparisons), ncol=3, byrow=TRUE))
      perm <- get_best_permutation(pair_comparisons_df)
      Q_list[[i]] <- Q_list[[i]][,perm]
      permutations[i,] <- perm
    }
  }

  #Check for bug in using the same column twice
  if (!all(unlist(lapply(Q_list, function(x){ length(unique(colnames(x)))==ncol(x) })))) stop("Duplicated column names in output Q matrices")

  #Rename columns
  column_names <- paste("Cluster ", seq(1,ncol(Q_list[[1]])))
  Q_list <- lapply(Q_list, function(x){
    colnames(x) <- column_names
    return(x)}
    )

  return(list(Q_list=Q_list, permutations=permutations))
}

G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}

get_best_permutation <- function(pair_df){
  #returns the best permutation of columns based on a dataframe of pairwise comparisons
  K <- max(pair_df[,2])
  best_pairs_df <- data.frame(x=rep(0.0, K), y=rep(0,K), z=rep(0,K))
  for (i in 1:K){
    best_pairs_df[i,] <- pair_df[which.max(pair_df[,1]),]
    dont_keep <- pair_df[,2]!=pair_df[which.max(pair_df[,1]),2]
    dont_keep <- dont_keep & (pair_df[,3]!=pair_df[which.max(pair_df[,1]),3])
    pair_df <- pair_df[dont_keep, ]
  }
  best_pairs_df <- best_pairs_df[order(best_pairs_df[,2]),]
  return(best_pairs_df[,3])
}

J_largeK <- function(Q_sub_list, Q_x, y, z){
  1/length(Q_sub_list) * sum(
    unlist(lapply(Q_sub_list, function(Q){
      G(Q[,y,drop=FALSE], Q_x[,z,drop=FALSE])
      }))
    )
}

J_perm <- function(perm, Q_x, Q_sub_list){
  J(Q_x[,perm], Q_sub_list)
}

J <- function(Q_x, Q_sub_list){
  sum(unlist(map_dbl(Q_sub_list, G, Q_x)))/length(Q_sub_list)
}

