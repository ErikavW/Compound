START_par <- function(N,K,D,M){
  #Startwaardes
  alpha<-matrix(1,nrow=K,ncol=D)

  beta<-matrix(NA,nrow=K,ncol=max(M))
  for (i in 1:K){
    beta[i,]<-sort(rnorm(max(M),0,runif(1,0.5,2)),decreasing=F)
  }
  beta<-cbind(-Inf,beta,Inf)  ## waardes voor 0 en M+1
  
  
  theta<-rmvnorm(N,c(0,0,0),matrix(c(1,0,0,0,1,0,0,0,1),nrow=D,ncol=D))
  
  return(list(alpha,theta,beta))
}
