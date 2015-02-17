CLEAN_itempar<-function(y,beta){
  yi<-y
  yi  <- apply(yi,2,function(x) { x[x==9] <- NA; x } )
  
  for (i in 1:K){
    if (any(!is.na(yi[,i]))){
      for (p in (1:M[i])[!(1:M[i] %in% yi[,i])]+1){
        beta[i,p]<-rtnorm(1,beta[i,p],sd.MH,left=beta[i,(p-1)],right=beta[i,(p+2)])
      }
    }else{
      beta[i,(1:max(M))+1]<-NA
    }
  }
}
