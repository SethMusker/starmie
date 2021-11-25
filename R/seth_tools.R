pairwiseSimilarityMatrix<-function(Q_list){
  R <- length(Q_list)
  combs<-combn(1:R,2)
  H <- matrix(nrow=R,ncol=R)
  for(i in 1:ncol(combs)){
    a<-combs[1,i]
    b<-combs[2,i]
    H[a,b]<-starmie:::G(Q_list[[a]],Q_list[[b]])
    H[b,a]<-H[a,b]
  }
  diag(H)<-1
  return(H)
}
environment(pairwiseSimilarityMatrix) <- asNamespace('starmie')
