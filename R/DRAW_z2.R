DRAW_z2 <- function(y,alpha,theta,beta,response,M){
  #install.packages("BayesBridge")
  library("BayesBridge") # Truncated normale verdleing
  
  N <- nrow(theta)
  D <- ncol(theta)
  K <- nrow(beta)
  
  z<-matrix(c(0),nrow=N,ncol=K)
  
  for (n in 1:N){
    for (i in 1:K){
      true<-sum(alpha[i,1:D]*theta[n,1:D]) # meer dimensies
      if (y[n,i]==9){
        # Impute missing data. Draw from normal distribution
        z[n,i]<-rnorm(1, mean=true, sd=1)
      }else{
        
        if (response == "polytoom"){
          if (y[n,i]==1){
            z[n,i]<- rtnorm(1,mu=true, right=beta[i,y[n,i]])
          }else{
            if (y[n,i]==M[i]){
              z[n,i]<- rtnorm(1,mu=true, left=beta[i,(y[n,i]-1)])
            }else{
              # Replace ordinal data with continues data. Trek uit normaal verdeling
              z[n,i]<- rtnorm(1,mu=true,left=beta[i,(y[n,i])-1], right=beta[i,y[n,i]])
            }
          }
        }else{
          true<-sum(alpha[i,1:D]*theta[n,1:D])-beta[i,1]
          if (y[n,i]==1){
            # Draw z from truncated normal left
            z[n,i]<-rtnorm(1, mu=true, sig=1, left=0)
          }else{
            if (y[n,i]==0){
              # Draw z from truncated normal rigth
              z[n,i]<-rtnorm(1, mu=true, sig=1, right=0)
            }else{
              print("Error DRAW_z")
            }
          }
        }
      }
    }
  }
  return(z)
}