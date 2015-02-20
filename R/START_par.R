START_par <- function(N,K,D,M){
  #Startwaardes
  alpha<-matrix(c(rep(1,K)),nrow=K,ncol=D, byrow=T)

  beta<-matrix(NA,nrow=K,ncol=max(M))
  for (i in 1:K){
    beta[i,]<-sort(rnorm(max(M),0,runif(1,0.5,2)),decreasing=F)
  }
  theta<-rmvnorm(N,rep(0,D),diag(D))
  
  return(list(alpha,theta,beta))
}
