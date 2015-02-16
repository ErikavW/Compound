#' Gibbs Simulator with data augmentation.
#' 
#' Estimates parameters of a GRM and the reliability and weigths of the most optimal compound score.
#' \deqn{P (X = j | \theta_1 ) = \Phi (\alpha \theta - \beta{ij}) - \Phi (\alpha \theta - \beta{i(j+1)})}
#' 
#' Details
#' 
#' @param y data matrix of \code{n x k}
#' @param alpha (optional) starting/fixed values for alpha otherwise draw.start should be true. 
#' @param theta (optional) starting/fixed values for theta otherwise draw.start should be true. 
#' @param beta (optional) starting/fixed values for beta otherwise draw.start should be true. 
#' @param response "dichotoom", "polytoom", "continue" 
#' @param Sim Amount of simulations
#' @param D Amount of dimensions
#' @param M Vector with the max amount of categories in each question.
#' @param identification 1 = Fikseren van eerste D items van alpha = a 0 0 a a 0 a a a,  Samen met sigma.theta = I, mu.theta = 0
#' @param draw.start True when starting values of 
#' @param schat.itemparameters True when estimating itemparameters is needed, False when itemparameters are given by alpha, theta and beta. 
#' 
#' @return Returns a list with the following estimatres:
#' @return $\code{alpha} factorloadings of the items.
#' @return \code{theta} person estimates of the items.
#' @return \code{beta} location of the items.
#' @return \code{rho} the reliability of the optimal scores 
#' @return \code{weights} the weigths for each of the dimensions
#' 
#' @export

compound<-function(y,alpha,theta,beta,response,Sim,D,M,identification=1,draw.start=T,fix.itemparameters=T){
  # Benodigde packages
  #install.packages("BayesBridge")
  library("BayesBridge") # Truncated normale verdleing
  library("mvtnorm")     # Multivariate normale verdeling
  library("MCMCpack")    # Inverse Wishart verdeling
  
  N<-nrow(y)
  K<-ncol(y)
  
  # beta <- CLEAN_itempar(y,beta)
  if (draw.start == T){
    source("C:\\Users\\Erika\\Documents\\R\\Optimal MCMC\\Data-archive\\Functions\\START_par.r")
    
    par <- START_par(N,K,D,M)
    alpha<-par[[1]]
    theta<-par[[2]]
    beta<-par[[3]]
  }
  
  # Prior waarden
  v0 <- c(-1) # Df inverse wishart
  k0 <- 0.01  # Weigth for prior theta verdeling
  mu0 <- 0    # Mean theta
  lam0 <- matrix(0.01,nrow=D,ncol=D) #Shape inverse wishart
  
  # Itemparameters priors
  mu.alpha<-rep(0,D+1)
  var.alpha<-rep(0.5,D+1)
  
  u.p<-rbind(c(mu.alpha))
  var.p<-diag(c(var.alpha))
  
  # Creeeren matrices die later gebruikt worden voor latente continue data
  rho<-matrix(NA,nrow=1,ncol=Sim)
  weights <-matrix(NA,nrow=D,ncol=Sim)
  
  # Creeeren matrices die later gebruikt worden voor creeren caterpillarplots alpha en beta
  A<-array(c(0),dim=c(K,D,Sim))
  B<-array(c(0),dim=c(K,max(M)+2,Sim))
  
  if (D==1){
    print("Only one dimension no weigths estimated, only person and itemparameters")
  }
  
  # Start iteraties
  for (s in 1:Sim){
    # Stap 1.Trekking Sigma en mu theta
    mu.sig <- DRAW_mu_sig(theta)
    
    # Stap 2. Trekken continue data voor y=0,1,..,M (9 is missing), uit truncated norm
    if (response == "continue" ){
      z <- (y-mean(y))/sd(y)
    }else{
      z <- DRAW_z(y,alpha,theta,beta,response)
    } 
    # Stap 3. Trek theta uit normaal verdeling mu en sigma theta.
    theta <- DRAW_theta(z=z,sig.theta=mu.sig[[2]],u=mu.sig[[1]],alpha=alpha,theta=theta)
    
    # Stap 4. Calculate optimal weights
    if (D !=1){
      rhox <- OPTIMAL_weights(true.A=mu.sig[[2]],obs.B=cov(theta),D=D)  
      rho[,s]<-rhox[[1]]
      weights[,s]<-rhox[[2]]      
    }
    
    
    # Stap 5 Trek itemparameters  
    if(fix.itemparameters==F){
      beta <- DRAW_beta(y,z,response,alpha,theta,beta,M,sd.MH=0.01, var.p,u.p)
      alpha <- DRAW_alpha(z,alpha,theta,var.p,u.p)
      
      # Matrixes met caterpillar voor alpha en beta
      A[,,s]<-alpha
      B[,,s]<-beta
    }
  }
  return(list(alpha=alpha,theta=theta,beta=beta[,2:max(M)+1],rho=rho,weights=weights))
}
