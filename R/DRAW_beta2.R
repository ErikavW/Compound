DRAW_beta2 <- function(y,z,response,alpha,theta,beta,M,sd.MH=0.05,var.p,u.p){
  b<-beta
  K<-nrow(beta)
  N<-nrow(y)
  
  if (response =="polytoom"){
    for (i in 1:K){
      temp1<-1
      temp2<-1
      
      for (j in 1:M[i]){
        if(j ==1){
           b[i,1]<-rtnorm(1,beta[i,1],sd.MH,right=beta[i,(2)]) #Trekking van b
        }else{
          if(j ==M[i]){
            b[i,M[i]]<-rtnorm(1,beta[i,M[i]],sd.MH,left=b[i,(M[i]-1)]) #Trekking van b
          }else{
            b[i,j]<-rtnorm(1,beta[i,j],sd.MH,left=b[i,(j-1)],right=beta[i,(j+1)]) #Trekking van b
          }
        }
      }
      
      for (n in 1:N){
        ee<- y[n,i]
        true <- alpha[i]*theta[n]
        if (y[n,i]!=0){
          if(y[n,i]==1){
            Ab<-1-pnorm(true-b[i,ee])
            Bb<-1-pnorm(true-beta[i,ee])
          }else{
            if(y[n,i]==M[i]){
              Ab<-pnorm(true-b[i,ee-1])
              Bb<-pnorm(true-beta[i,ee-1])
            }else{
              Ab<-pnorm(true-b[i,ee-1])-pnorm(true-b[i,ee])
              Bb<-pnorm(true-beta[i,ee-1])-pnorm(true-beta[i,ee])
            }          
          }
        }
        temp1<-temp1*(Ab/Bb)   
      }
      for (n in 1:N){
        for(j in 1:M[i]){
          if(j==1){
            Ab<-pnorm((beta[i,2]-beta[i,1]),sd.MH)
            Bb<-pnorm((b[i,2]-b[i,1]),sd.MH)
          }else{
            if(j==M[i]){
              Ab<-1-pnorm((b[i,(M[i]-1)]-beta[i,M[i]]),sd.MH)
              Bb<-1-pnorm((beta[i,(M[i]-1)]-b[i,M[i]]),sd.MH)
            }else{
              Ab<-pnorm((beta[i,(j+1)]-beta[i,j]),sd.MH)-pnorm((b[i,(j-1)]-beta[i,j]),sd.MH)
              Bb<-pnorm((b[i,(j+1)]-b[i,j]),sd.MH)-pnorm((beta[i,(j-1)]-b[i,j]),sd.MH)
            }
          }
        }
        temp2<-temp2*(Ab/Bb)        
      } 
      R=temp1*temp2
      ru<-runif(1,0,1)
      if (R<=ru){
        beta[i,1:M[i]]<-b[i,M[i]]
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
      beta[i,1]<-YY[1] 
    }
  }  
  return(beta)
}
