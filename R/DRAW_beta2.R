DRAW_beta <- function(y,z,response,alpha,theta,beta,M,sd.MH=0.05, var.p,u.p){
  b<-beta
  K<-nrow(beta)
  N<-nrow(y)
  
  if (response =="polytoom"){
    for (i in 1:K){
      temp1<-1
      temp2<-1
      
      if(j ==1){
        b[i,j]<-rtnorm(1,beta[i,j],sd.MH,right=beta[i,(j+1)]) #Trekking van b
      }
      if(j ==M[i]){
        b[i,j]<-rtnorm(1,beta[i,j],sd.MH,left=b[i,(j-1)]) #Trekking van b
      }
        
      for (j in 2:(M[i]-1)){
        b[i,j]<-rtnorm(1,beta[i,j],sd.MH,left=b[i,(j-1)],right=beta[i,(j+1)]) #Trekking van b
      }
      
      for (n in 1:N){
        ee<- y[n,i]
        true <- alpha[i]*theta[n]
        if(y[n,i]=1){
          Ab<-1-pnorm(true-b[i,ee])
          Bb<-1-pnorm(true-beta[i,ee])
        }else{
          if(y[n,i]=max(m)){
            Ab<-pnorm(true-b[i,ee-1])
            Bb<-pnorm(true-beta[i,ee-1])
          }else{
            Ab<-pnorm(true-b[i,ee-1])-pnorm(true-b[i,ee])
            Bb<-pnorm(true-beta[i,ee-1])-pnorm(true-beta[i,ee])
          }          
        }
        temp1<-temp1*(Ab/Bb)        
      }
      for (n in 1:N){
        if(j==1){
          Ab<-pnorm((beta[i,(j+1)]-beta[i,j]),sd.MH)-1
          Bb<-pnorm((b[i,(j+1)]-b[i,j]),sd.MH)-1
        }
        if(j==M[i]){
          Ab<-0-pnorm((b[i,(j-1)]-beta[i,j]),sd.MH)
          Bb<-0-pnorm((beta[i,(j-1)]-b[i,j]),sd.MH)
        }
        for (j in 2:(M[i]-1)){
          Ab<-pnorm((beta[i,(j+1)]-beta[i,j]),sd.MH)-pnorm((b[i,(j-1)]-beta[i,j]),sd.MH)
          Bb<-pnorm((b[i,(j+1)]-b[i,j]),sd.MH)-pnorm((beta[i,(j-1)]-b[i,j]),sd.MH)
        }          
      temp2<-temp2*(Ab/Bb)        
       
      } 
      R=temp1*temp2
      ru<-runif(1,0,1)
      if (R<=ru){
        beta[i,1:M[i])]<-b[i,M[i]]
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
