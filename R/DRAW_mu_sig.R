DRAW_mu_sig<-function(theta, identification = 1, v0=v0, k0=k0, mu0=mu0,lam0=lam0){
  
  if (identification == 1){
    # Mean theta is 0 <-- Identificatie
    u<-matrix(c(0),ncol=ncol(theta))
    
    # Sigma theta is identiteitsmatrix <-- Identificatie
    sig.theta <- diag(ncol(theta))
    
    
  }else{
    N <- nrow(theta)
    kn <- k0+N
    vn <- v0+N
    lam <- matrix(c(lam0),nrow=ncol(theta),ncol=ncol(theta))
    
    # Trekking mean theta
    theta.mean <- apply(theta,2,mean)
    mun <- k0/(k0+n)*mu0+n/(k0+n)*theta.mean
    u <- rmvnorm(1,mean=mun,sigma=(sig.theta/kn))
    
    # Trekking Sigma theta
    S<-cov(theta)*(N-1)   #Sum of squares
    lamn<-lam+S+(k0*n)/(k0+n)*(theta.mean-mu0)%*%t(theta.mean-mu0) #Shape for inverse wishart
    sig.theta<-riwish(vn,lamn)  #trekking Sigma theta uit inverse wishart verdeling (niet solve(lamn), zoals in artikel beguin)   
    
  }
  return(list(u,sig.theta))
}
  




