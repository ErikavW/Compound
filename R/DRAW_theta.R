DRAW_theta<-function(z,sig.theta,u,alpha,theta){
  D <- ncol(theta)
  N <- nrow(theta)
  # Creeer "sd" van theta.
  L <- chol(sig.theta)
  B <- alpha%*%L
  BB <- t(B)%*%B
  sig <- solve(BB)
  theta.u <- solve(diag(D)+BB) # theta.u is de varcovar van de orthogonally standardized ability.

  for (n in 1:N){ 
    theta.0s<-sig%*%t(B)%*%(z[n,]-alpha%*%t(u)) # Beta zit in z nog verwerkt
    theta.0<-rmvnorm(1,(theta.u%*%BB%*%theta.0s),theta.u)
    theta[n,]<-L%*%t(theta.0)+t(u)
  }
  return(theta)
}

