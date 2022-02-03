SIMEXBoost<-function(Y,Xstar,zeta=c(0,0.25,0.5,0.75,1),B=500,type="normal",sigmae,Iter=100,Lambda=0,Extrapolation="linear"){

  if(length(zeta)<1) return("ERROR:zeta length should larger than 1")
  if(B<1) return("ERROR:B should larger than 1")
  if(dim(Xstar)[2]!=dim(sigmae)[2]) return("ERROR:sigmae should have the same dimention with Xstar")
  if(Lambda<(-0.00001)) return("ERROR:Lambda should larger than 0")
  if(type!="normal"&type!="binary"&type!="poisson"&type!="AFT-normal"&type!="AFT-loggamma") return("ERROR:tpye set error")
  p=dim(Xstar)[2];n=dim(Xstar)[1];betazero=t(rep(0,p));
  if(type=="normal"){
    betarbind<-NULL;
    for(cc in 1:length(zeta)){
      betaihat<-NULL
      for(i in 1:B){
        ec<-mvrnorm(n=n,rep(0,p),sigmae)
        Wetabt<-Xstar+sqrt(zeta[cc])*ec#continous
        Wetab<-t(Wetabt);time=0;
        for(j in 1:Iter){
          betazero1=betazero
          ubetastep<-NULL
          for(i in 1:length(Y)){
            ubeta<-(Y[i]-sum(Wetab[,i]*betazero1))
            ubeta1<-Wetab[,i]*ubeta - Lambda * betazero1
            ubetastep<-rbind(ubetastep,ubeta1)
          }
          #¼ÆËãubeta
          u<-NULL
          for(i in 1:p){
            u<-c(u,sum(ubetastep[,i])/length(Y))
          }
          for(i in 1:p){
            Delta=u
            #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
            if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
            #if(time==Iter-1) betazero[which(abs(betazero)<0.16)]=0
          }
          time=time+1
        }
        betaihat<-rbind(betaihat,betazero)
      }
      betarbind<-rbind(betarbind,betaihat)
    }
    i=0;meanBi<-NULL
    while((i*B)!=dim(betarbind)[1]){
      Bi<-colMeans(betarbind[((i*B)+1):((i+1)*B),]);
      meanBi<-rbind(meanBi,Bi)
      i=i+1
    }
    for(i in 1:(length(zeta)*p)){
      if(abs(meanBi[i])<0.21) meanBi[i]=0
    }
    if(Extrapolation=="linear"){
      gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x)
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        #gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        #gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat#+gamma2hat
    }
    if(Extrapolation=="quadratic"){
      gamma1hat<-NULL;gamma0hat<-NULL;gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x+I(x*x))
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat+gamma2hat
    }
    for(i in 1:p){
      if(abs(betacorrect[i])<0.1) betacorrect[i]=0
    }
  }
  if(type=="binary"){
    betarbind<-NULL;
    for(cc in 1:length(zeta)){
      betaihat<-NULL
      for(i in 1:B){
        ec<-mvrnorm(n=n,rep(0,p),sigmae)
        Wetabt<-Xstar+sqrt(zeta[cc])*ec
        Wetab<-t(Wetabt);time=0;
        for(j in 1:Iter){
          betazero1=betazero
          ubetastep<-NULL
          for(i in 1:length(Y)){
            ubeta<-Y[i]-(exp(sum(Wetab[,i]*betazero1))/(1+exp(sum(Wetab[,i]*betazero1))))
            ubeta1<-Wetab[,i]*ubeta - Lambda * betazero1
            ubetastep<-rbind(ubetastep,ubeta1)
          }
          #¼ÆËãubeta
          u<-NULL
          for(i in 1:p){
            u<-c(u,sum(ubetastep[,i])/length(Y))
          }
          for(i in 1:p){
            Delta=u
            #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
            if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
          }
          time=time+1
        }
        betaihat<-rbind(betaihat,betazero)
      }
      betarbind<-rbind(betarbind,betaihat)
    }
    i=0;meanBi<-NULL
    while((i*B)!=dim(betarbind)[1]){
      Bi<-colMeans(betarbind[((i*B)+1):((i+1)*B),]);
      meanBi<-rbind(meanBi,Bi)
      i=i+1
    }
    for(i in 1:(length(zeta)*p)){
      if(abs(meanBi[i])<0.21) meanBi[i]=0
    }
    if(Extrapolation=="linear"){
      gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x)
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        #gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        #gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat#+gamma2hat
    }
    if(Extrapolation=="quadratic"){
      gamma1hat<-NULL;gamma0hat<-NULL;gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x+I(x*x))
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat+gamma2hat
    }
    for(i in 1:p){
      if(abs(betacorrect[i])<0.21) betacorrect[i]=0
    }
  }
  if(type=="poisson"){
    betarbind<-NULL;
    for(cc in 1:length(zeta)){
      betaihat<-NULL
      for(i in 1:B){
        ec<-mvrnorm(n=n,rep(0,p),sigmae)
        Wetabt<-Xstar+sqrt(zeta[cc])*ec
        Wetab<-t(Wetabt);time=0;
        for(j in 1:Iter){
          betazero1=betazero
          ubetastep<-NULL
          for(i in 1:length(Y)){
            ubeta<-(Y[i]-exp(sum(Wetab[,i]*betazero1)))
            ubeta1<-Wetab[,i]*ubeta - Lambda * betazero1
            ubetastep<-rbind(ubetastep,ubeta1)
          }
          #¼ÆËãubeta
          u<-NULL
          for(i in 1:p){
            u<-c(u,sum(ubetastep[,i])/length(Y))
          }
          for(i in 1:p){
            Delta=u
            #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
            if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
          }
          time=time+1
        }
        betaihat<-rbind(betaihat,betazero)
      }
      betarbind<-rbind(betarbind,betaihat)
    }
    i=0;meanBi<-NULL
    while((i*B)!=dim(betarbind)[1]){
      Bi<-colMeans(betarbind[((i*B)+1):((i+1)*B),]);
      meanBi<-rbind(meanBi,Bi)
      i=i+1
    }
    for(i in 1:(length(zeta)*p)){
      if(abs(meanBi[i])<0.21) meanBi[i]=0
    }
    if(Extrapolation=="linear"){
      gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x)
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        #gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        #gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat#+gamma2hat
    }
    if(Extrapolation=="quadratic"){
      gamma1hat<-NULL;gamma0hat<-NULL;gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x+I(x*x))
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat+gamma2hat
    }
    for(i in 1:p){
      if(abs(betacorrect[i])<0.1) betacorrect[i]=0
    }
  }
  if(type=="AFT-normal"){
    betarbind<-NULL;YL=Y[,1];YR=Y[,2]
    for(cc in 1:length(zeta)){
      betaihat<-NULL;
      TT1=YL[which(YL!=0)];TT2=YR[which(YL!=0)]
      interval1=which(YL==0)
      interval=which(TT1==(YR[-interval1]))
      TT1real<-TT1[interval]
      TT1censor<-TT1[-interval]
      R<-TT2[-interval]
      LR<-cbind(TT1censor,R)

      for(i in 1:B){
        ec<-mvrnorm(n=n,rep(0,p),sigmae)
        Wetabt<-Xstar+sqrt(zeta[cc])*ec
        Wetab<-t(Wetabt);time=0;
        X1<-Wetab[,-interval1]
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
          for(i in 1:p){
            Delta=u
            #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
            if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
            if(time==Iter-1) betazero[which(abs(betazero)<0.41)]=0
          }
          time=time+1
        }
        betaihat<-rbind(betaihat,betazero)
      }
      betarbind<-rbind(betarbind,betaihat)
    }
    i=0;meanBi<-NULL
    while((i*B)!=dim(betarbind)[1]){
      Bi<-colMeans(betarbind[((i*B)+1):((i+1)*B),]);
      meanBi<-rbind(meanBi,Bi)
      i=i+1
    }
    for(i in 1:(length(zeta)*p)){
      if(abs(meanBi[i])<0.41) meanBi[i]=0
    }
    if(Extrapolation=="linear"){
      gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x)
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        #gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        #gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat#+gamma2hat
    }
    if(Extrapolation=="quadratic"){
      gamma1hat<-NULL;gamma0hat<-NULL;gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x+I(x*x))
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat+gamma2hat
    }
    for(i in 1:p){
      if(abs(betacorrect[i])<0.41) betacorrect[i]=0
    }
  }
  if(type=="AFT-loggamma"){
    betarbind<-NULL;YL=Y[,1];YR=Y[,2]
    for(cc in 1:length(zeta)){
      betaihat<-NULL;
      TT1=YL[which(YL!=0)];TT2=YR[which(YL!=0)]
      interval1=which(YL==0)
      interval=which(TT1==(YR[-interval1]))
      TT1real<-TT1[interval]
      TT1censor<-TT1[-interval]
      R<-TT2[-interval]
      LR<-cbind(TT1censor,R)

      for(i in 1:B){
        ec<-mvrnorm(n=n,rep(0,p),sigmae)
        Wetabt<-Xstar+sqrt(zeta[cc])*ec
        Wetab<-t(Wetabt);time=0;
        X1<-Wetab[,-interval1]
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
          for(i in 1:p){
            Delta=u
            #if(abs(Delta[i])>0.6*max(abs(Delta))) betazero[i]=betazero[i]+0.01*sign(Delta[i])#Delta[i]=Delta[i]+0.05*sign(Delta[i])
            if(abs(Delta[i])>0.9*max(abs(Delta))) betazero[i]=betazero[i]+0.05*sign(Delta[i])
            if(time==Iter-1) betazero[which(abs(betazero)<0.41)]=0
          }
          time=time+1
        }
        betaihat<-rbind(betaihat,betazero)
      }
      betarbind<-rbind(betarbind,betaihat)
    }
    i=0;meanBi<-NULL
    while((i*B)!=dim(betarbind)[1]){
      Bi<-colMeans(betarbind[((i*B)+1):((i+1)*B),]);
      meanBi<-rbind(meanBi,Bi)
      i=i+1
    }
    for(i in 1:(length(zeta)*p)){
      if(abs(meanBi[i])<0.41) meanBi[i]=0
    }
    if(Extrapolation=="linear"){
      gamma1hat<-NULL;gamma0hat<-NULL;#gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x)
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        #gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        #gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat#+gamma2hat
    }
    if(Extrapolation=="quadratic"){
      gamma1hat<-NULL;gamma0hat<-NULL;gamma2hat<-NULL;
      for(i in 1:p){
        x<-zeta
        y<-NULL
        for(j in 1:length(zeta)){
          y1<-meanBi[j,i]
          y<-c(y,y1)
        }
        c<-lm(y~x+I(x*x))
        gamma0<-as.numeric(c$coefficients[1])
        gamma1<-as.numeric(c$coefficients[2])
        gamma2<-as.numeric(c$coefficients[3])
        gamma0hat<-c(gamma0hat,gamma0)
        gamma1hat<-c(gamma1hat,gamma1)
        gamma2hat<-c(gamma2hat,gamma2)
      }
      betacorrect<-gamma0hat-gamma1hat+gamma2hat
    }
    for(i in 1:p){
      if(abs(betacorrect[i])<0.41) betacorrect[i]=0
    }

  }
  return(list("BetaHatCorrect"=betacorrect))
}

