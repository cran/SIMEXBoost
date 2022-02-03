ME_Data<-function(X,beta,type="normal",sigmae,pr0=0.5){
  if(dim(X)[2]!=dim(sigmae)[2]) return("ERROR:sigmae should have the same dimention with Xstar")
  if(length(beta)!=dim(X)[2]) return("ERROR:length(beta) should equal to X de dimension")
  if(type!="normal"&type!="binary"&type!="poisson"&type!="AFT-normal"&type!="AFT-loggamma") return("ERROR:tpye set error")
  if(pr0<0 | pr0>1) return("ERROR:pr0 should in (0,1)")
  p=dim(X)[2];n=dim(X)[1]
  xtbeta<-(t(X))*beta
  Y=NULL
  #compute Y
  if(type=="normal"){
    for(i in 1:n){
      e<-rnorm(1)
      y=sum(xtbeta[,i])+e
      Y=cbind(Y,y)
    }
  }
  if(type=="binary"){
    for(i in 1:n){
      y=rbinom(1,1,exp(sum(xtbeta[,i]))/(1+exp(sum(xtbeta[,i]))))
      Y=cbind(Y,y)
    }
  }
  if(type=="poisson"){
    for(i in 1:n){
      e<-rnorm(1)
      y=rpois(1,exp(sum(xtbeta[,i])))
      Y=cbind(Y,y)
    }
  }
  if(type=="AFT-normal"){
    logT=NULL
    for(i in 1:n){
      e<-rnorm(1)
      logt=sum(xtbeta[,i])+e
      logT=cbind(logT,logt)
    }
    Y=exp(logT)
    A<-runif(length(Y),0,sort(Y)[0.95*n])
    for(i in 1:n){
      if(A[i]>Y[i]) Y[i]=0
    }
    Y1=Y[Y!=0]
    delta=rbinom(length(Y1),1,pr0)
    interval=which(delta==1)
    Y1real<-Y1[interval]
    Y1censor<-Y1[-interval]
    interval2=which(Y==0)
    A1<-A[-interval2]
    A1real<-A1[interval]
    A1censor<-A1[-interval]

    LR=NULL
    for(i in 1:length(Y1censor)){
      u=A1censor[i]
      u1=NULL
      while(u<=Y1censor[i]){
        u=u+0.1+runif(1,0,1)
        u1=c(u1,u)
      }
      LR<-rbind(LR,c(max(A1censor[i],u1[(length(u1)-1)]),u1[length(u1)]))
    }
    YL<-Y;YR<-Y
    for(i in 1:length(Y1censor)){
      YL[which(Y1censor[i]==Y)]=LR[i,1]
      YR[which(Y1censor[i]==Y)]=LR[i,2]
    }
    Y<-rbind(YL,YR)
  }
  if(type=="AFT-loggamma"){
    logT=NULL
    for(i in 1:n){
      e<-log(rexp(1))+0.5772156649
      logt=sum(xtbeta[,i])+e
      logT=cbind(logT,logt)
    }
    Y=exp(logT)
    A<-runif(length(Y),0,sort(Y)[0.95*n])
    for(i in 1:n){
      if(A[i]>Y[i]) Y[i]=0
    }
    Y1=Y[Y!=0]
    delta=rbinom(length(Y1),1,pr0)
    interval=which(delta==1)
    Y1real<-Y1[interval]
    Y1censor<-Y1[-interval]
    interval2=which(Y==0)
    A1<-A[-interval2]
    A1real<-A1[interval]
    A1censor<-A1[-interval]

    LR=NULL
    for(i in 1:length(Y1censor)){
      u=A1censor[i]
      u1=NULL
      while(u<=Y1censor[i]){
        u=u+0.1+runif(1,0,1)
        u1=c(u1,u)
      }
      LR<-rbind(LR,c(max(A1censor[i],u1[(length(u1)-1)]),u1[length(u1)]))
    }
    YL<-Y;YR<-Y
    for(i in 1:length(Y1censor)){
      YL[which(Y1censor[i]==Y)]=LR[i,1]
      YR[which(Y1censor[i]==Y)]=LR[i,2]
    }
    Y<-rbind(YL,YR)
  }
  Xstar<-X+mvrnorm(n=n,rep(0,p),sigmae)
  List = list("response" = t(Y), "ME_covariate" = Xstar)
  return(List)
}
