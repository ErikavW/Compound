DRAW_alpha<-function(z,alpha,theta,var.p,u.p, identification = 1){
  N<-nrow(z)
  K<-ncol(z)
  D<-ncol(theta)
  X<-cbind(-1,theta) # Creeer vector regressiecoefficient
  XX<-solve(solve(var.p)+t(X)%*%X)      # 
  
  if (identification == 1){
    ###################### Manier 1 #########################
    # Fikseren van eerste D items van alpha = a 0 0         #
    #                                         a a 0         #
    #                                         a a a         # 
    # Samen met sigma.theta = I                             #
    # mu.theta = 0                                          #
    #########################################################
    
    for (i in 1:K){
      zeta<-XX%*%(solve(var.p)%*%t(u.p)+t(X)%*%z[,i])
      YY<-c(rmvnorm(1, mean=zeta,sigma=XX))
      while ((YY[1]|YY[2])<0){
        YY<-c(rmvnorm(1, mean=zeta,sigma=XX))
        tel<-tel+1
        if (tel>200){
          print(c(n,i))
        }
      }
      alpha[i,1:min(i,D)]<-YY[2:(min(i,D)+1)]
      alpha[max(i,D+1),]<-YY[2:(D+1)]
    }
  }else{
    ###################### Manier 2 #########################
    # Fikseren van eerste D items van alpha = I --> diag(D).# 
    # Sigma.theta is vrij                                   #
    # mu.theta = 0.                                         #
    #########################################################
    for (i in 1:K){
      zeta<-XX%*%(solve(var.p)%*%t(u.p)+t(X)%*%z[,i])    # mean van itemparameters
      YY<-c(rmvnorm(1, mean=zeta,sigma=XX))
      while ((YY[1]|YY[2])<0){
        YY<-c(rmvnorm(1, mean=zeta,sigma=XX))
        tel<-tel+1
        if (tel>200){
          print(c(n,i))
        }
      }
      alpha[max(i,D+1),]<-YY[2:D+1]   # Alleen wegschrijven alpha na eerste D items, eerste D --> I.
    } 
  }
  return(alpha)
}

