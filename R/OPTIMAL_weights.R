OPTIMAL_weights <- function(true.A,obs.B,D){
  # Optimal Weigths 
  delta <- diag(eigen(true.A)$values)
  U1 <- eigen(true.A)$vectors                   #Non-zero eigenvalues
  
  C11 <- solve(sqrt(delta)) %*% t(U1) %*% obs.B %*% U1 %*% solve(sqrt(delta))
  
  Gamma <- diag(eigen(C11)$values)
  V <- eigen(C11)$vectors
  
  Trans <- U1 %*% solve(sqrt(delta)) %*% V
  
  #zx<-sqrt(1/Gamma[D,D])
  zx<-matrix(c(sqrt(1/Gamma[1,1]),rep(0,(D-1))),nrow=D)
  x<- (Trans) %*% zx
  
  rho_e<-(t(x)%*% true.A %*% x) /(t(x)%*% obs.B %*% x) #Betrouwbaarheid zo groot mogelijk              
  
  rho<-(t(x)%*% true.A %*% x)/(t(x)%*% obs.B %*% x)
  weights<-x
  ## End Optimal Weigths
  
  return(list(rho,x))
}
  