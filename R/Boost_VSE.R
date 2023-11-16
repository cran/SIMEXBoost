Boost_VSE<-function(Y,Xstar,type="normal",Iter=200,Lambda=0){
  if(type!="normal"&type!="binary"&type!="poisson"&type!="AFT-normal"&type!="AFT-loggamma") return("ERROR:tpye set error")
  p=dim(Xstar)[2];n=dim(Xstar)[1];betazero=rep(0,p);X=Xstar;time=0
  if(type=="normal"){
    for(j in 1:Iter){
      betazero1=betazero
      ubetastep<-t(X) %*% (Y - X %*% betazero1) - Lambda * betazero1


      Delta<-ubetastep / length(Y)
      Jset = which(abs(Delta)>0.9*max(abs(Delta)))
      betazero[Jset]=betazero[Jset]+0.05*sign(Delta[Jset])

      #for(i in 1:p){
      #  Delta=u
      #  if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
      #}
      time=time+1
    }
  }
  if(type=="binary"){
    for(j in 1:Iter){
      betazero1=betazero
      ubetastep<-t(X) %*% (Y - (exp(X %*% betazero1)/(1+exp(X %*% betazero1)))   ) - Lambda * betazero1

      Delta<-ubetastep / length(Y)
      Jset = which(abs(Delta)>0.9*max(abs(Delta)))
      betazero[Jset]=betazero[Jset]+0.05*sign(Delta[Jset])
      #for(i in 1:p){
      #  Delta=u
      #  if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
      #}
      time=time+1
    }


  }
  if(type=="poisson"){
    for(j in 1:Iter){
      betazero1=betazero
      ubetastep<-t(X) %*% (Y - exp(X %*% betazero1)) - Lambda * betazero1


      Delta<-ubetastep / length(Y)
      Jset = which(abs(Delta)>0.9*max(abs(Delta)))
      betazero[Jset]=betazero[Jset]+0.05*sign(Delta[Jset])
      #for(i in 1:p){
      #  Delta=u
      #  if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
      #}
      time=time+1
    }


  }
  if(type=="AFT-normal"){
    YL=Y[,1];YR=Y[,2];
    TT1=YL[which(YL!=0)];TT2=YR[which(YL!=0)]
    interval1=which(YL==0)
    interval=which(TT1==(YR[-interval1]))
    TT1real<-TT1[interval]
    TT1censor<-TT1[-interval]
    R<-TT2[-interval]
    LR<-cbind(TT1censor,R)
    X1<-X[,-interval1]
    X1real<-X1[,interval]
    X1censor<-X1[,-interval]
    logTT1real<-log(TT1real);logTT1censor<-log(TT1censor)
    functionx<-function(x){
      b=1/(exp(1)*(exp(1)-1))
      a=(log(x)/x)*(b/((pi)*(((log(x)/x)^2)+b^2)))
      return(a)
    }
    function2<-function(x){
      b=1/(exp(1)*(exp(1)-1))
      a=(b/((pi)*(((log(x)/x)^2)+b^2)))
      return(a)
    }
    for(j in 1:Iter){
      betazero1=betazero
      Lbeta1<-NULL
      for(i in 1:length(TT1censor)){
        Lbeta<-LR[,1][i]*exp(-sum(X1censor[,i]*betazero1))
        Lbeta1<-cbind(Lbeta1,Lbeta)
      }
      Rbeta1<-NULL
      for(i in 1:length(TT1censor)){
        Rbeta<-LR[,2][i]*exp(-sum(X1censor[,i]*betazero1))
        Rbeta1<-cbind(Rbeta1,Rbeta)
      }
      ubetastep1<-NULL
      for(i in 1:length(logTT1real)){
        ubeta<-(logTT1real[i]-sum(X1real[,i]*betazero1))/
          (TT1real[i]*exp(-sum(X1real[,i]*betazero1)))
        ubeta1<-X1real[,i]*ubeta-Lambda*betazero1
        ubetastep1<-rbind(ubetastep1,ubeta1)
      }
      dFtgeneration<-NULL;
      for(i in 1:length(TT1censor)){
        a<-runif(5000,Lbeta1[i],Rbeta1[i])
        a1<-functionx(a)
        dFt=sum((Rbeta1[i]-Lbeta1[i])*a1)/length(a1)
        b<-runif(5000,Lbeta1[i],Rbeta1[i])
        b1<-function2(b)
        dFt1<-sum((Rbeta1[i]-Lbeta1[i])*b1)/length(b1)
        censorYi<-dFt/dFt1
        dFtgeneration<-rbind(dFtgeneration,censorYi)
      }
      for(i in 1:length(dFtgeneration)){
        if(dFtgeneration[i]=="Inf") dFtgeneration[i]=0
        if(dFtgeneration[i]=="-Inf") dFtgeneration[i]=0
        if(dFtgeneration[i]=="NaN") dFtgeneration[i]=0
      }
      ubetastep2<-NULL
      for(i in 1:length(dFtgeneration)){
        a<-X1censor[,i]*dFtgeneration[i]-Lambda*betazero1
        ubetastep2<-rbind(ubetastep2,a)
      }
      #compute u
      ubetastep<-rbind(ubetastep1,ubetastep2)
      u<-NULL
      for(i in 1:p){
        u<-c(u,sum(ubetastep[,i])/length(TT1))
      }
      Delta=u
      Jset = which(abs(Delta)>0.9*max(abs(Delta)))
      betazero[Jset]=betazero[Jset]+0.05*sign(Delta[Jset])
      if(time==Iter-1) betazero[which(abs(betazero)<0.11)]=0
      #for(i in 1:p){

      #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
      # if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])

      #}
      time=time+1
    }

  }
  if(type=="AFT-loggamma"){
    YL=Y[,1];YR=Y[,2];
    TT1=YL[which(YL!=0)];TT2=YR[which(YL!=0)]
    interval1=which(YL==0)
    interval=which(TT1==(YR[-interval1]))
    TT1real<-TT1[interval]
    TT1censor<-TT1[-interval]
    R<-TT2[-interval]
    LR<-cbind(TT1censor,R)
    X1<-X[,-interval1]
    X1real<-X1[,interval]
    X1censor<-X1[,-interval]
    logTT1real<-log(TT1real);logTT1censor<-log(TT1censor)
    functionx<-function(x){
      a=(log(x)/(x^3))
      return(a)
    }
    function2<-function(x){
      a=1/(x^2)
      return(a)
    }
    for(j in 1:Iter){
      betazero1=betazero
      Lbeta1<-NULL
      for(i in 1:length(TT1censor)){
        Lbeta<-LR[,1][i]*exp(-sum(X1censor[,i]*betazero1))
        Lbeta1<-cbind(Lbeta1,Lbeta)
      }
      Rbeta1<-NULL
      for(i in 1:length(TT1censor)){
        Rbeta<-LR[,2][i]*exp(-sum(X1censor[,i]*betazero1))
        Rbeta1<-cbind(Rbeta1,Rbeta)
      }
      ubetastep1<-NULL
      for(i in 1:length(logTT1real)){
        ubeta<-(logTT1real[i]-sum(X1real[,i]*betazero1))/
          (TT1real[i]*exp(-sum(X1real[,i]*betazero1)))
        ubeta1<-X1real[,i]*ubeta-Lambda*betazero1
        ubetastep1<-rbind(ubetastep1,ubeta1)
      }
      dFtgeneration<-NULL;
      for(i in 1:length(TT1censor)){
        a<-runif(5000,Lbeta1[i],Rbeta1[i])
        a1<-functionx(a)
        dFt=sum((Rbeta1[i]-Lbeta1[i])*a1)/length(a1)
        b<-runif(5000,Lbeta1[i],Rbeta1[i])
        b1<-function2(b)
        dFt1<-sum((Rbeta1[i]-Lbeta1[i])*b1)/length(b1)
        censorYi<-dFt/dFt1
        dFtgeneration<-rbind(dFtgeneration,censorYi)
      }
      for(i in 1:length(dFtgeneration)){
        if(dFtgeneration[i]=="Inf") dFtgeneration[i]=0
        if(dFtgeneration[i]=="-Inf") dFtgeneration[i]=0
        if(dFtgeneration[i]=="NaN") dFtgeneration[i]=0
      }
      ubetastep2<-NULL
      for(i in 1:length(dFtgeneration)){
        a<-X1censor[,i]*dFtgeneration[i]-Lambda*betazero1
        ubetastep2<-rbind(ubetastep2,a)
      }
      #compute u
      ubetastep<-rbind(ubetastep1,ubetastep2)
      u<-NULL
      for(i in 1:p){
        u<-c(u,sum(ubetastep[,i])/length(TT1))
      }
      Delta=u
      Jset = which(abs(Delta)>0.9*max(abs(Delta)))
      betazero[Jset]=betazero[Jset]+0.05*sign(Delta[Jset])
      if(time==Iter-1) betazero[which(abs(betazero)<0.11)]=0
      #for(i in 1:p){

      #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
      # if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
      # if(time==Iter-1) betazero[which(abs(betazero)<0.11)]=0
      #}
      time=time+1
    }


  }
  return(list("BetaHat"=betazero))
}

