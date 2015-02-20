DRAW_beta <- function(y,z,response,alpha,theta,beta,M,sd.MH=0.05, var.p,u.p){
  b<-beta
  K<-nrow(beta)
  N<-nrow(y)
  
  if (response =="polytoom"){
    for (i in 1:K){
      temp1<-1
      temp2<-1
      
      for (j in 2:M[i]){  #+1??
        b[i,j]<-rtnorm(1,beta[i,j],sd.MH,left=b[i,(j-1)],right=beta[i,(j+1)]) #Trekking van b
      }
      
      # Calculate acceptence rate, according to Johnson/Alber and the translation to irt of Glas in Ostini.
      for (j in 3:M[i]){  #+1??
        mean<-beta(i,j)
        Ab<-pnorm((mean-beta[i,(j+1)]),sd.MH)-pnorm((mean-b[i,(j-1)]),sd.MH) # Normaal Ogive model 
        Bb<-pnorm((mean-b[i,(j+1)]),sd.MH)-pnorm((mean-beta[i,(j-1)]),sd.MH)
        if (Bb==0 & Ab==0){  
          print("Bb & Ab are 0") 
          temp1<-temp1*1
        }else{
          temp1<-temp1*(Ab/Bb)
        }
      }
      for (n in 1:N){
        mean<-sum(alpha[i,]*theta[n,])
        
        if(y[n,i]==0){
          
        }else{
          ee<- y[n,i]+1
          Ab<-pnorm(mean-b[i,ee])-pnorm(mean-b[i,(ee-1)])
          Bb<-pnorm(mean-beta[i,ee])-pnorm(mean-beta[i,(ee-1)])
          if (Bb==0 & Ab==0){                         ### Wat als Bb en Ab beide 0 zijn???
            print("Bb & Ab are 0") 
            temp2<-temp2*1
          }else{
            temp2<-temp2*(Ab/Bb)
          }
        }
      }
      R=temp1*temp2
      ru<-runif(1,0,1)
      if (R<=ru){
        beta[i,2:(M[i]+1)]<-b[i,2:(M[i]+1)]
      }
    }    
  }else{
    X<-cbind(-1,theta)                    # Creeer vector regressiecoefficient
    XX<-solve(solve(var.p)+t(X)%*%X)      #  
    
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
      beta[i,2]<-YY[1] 
    }
  }
  
  return(beta)
}
