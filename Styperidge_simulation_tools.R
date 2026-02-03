
#Ege University Department of Statistics

library(MASS)
library("MASS", lib.loc="C:/Program Files/R/R-3.4.3/library")

meanbiasvarmse=function(est,parameter) {

  if (is.vector(est)){
    n=length(est)

    meanest=mean(est)
    biasest=meanest-parameter
    varest=var(est)
    mseest=biasest^2+varest
    y=list(meanest=meanest,biasest=biasest,varest=varest,mseest=mseest)
    return(y)
  }

  dimest=dim(est)
  c=dimest[2]

  if(c==1){
    meanest=colMeans(est)
    biasest=meanest-parameter
    varest=diag(var(est))
    mseest=biasest^2+varest
    y=list(meanest=meanest,biasest=biasest,varest=varest,mseest=mseest)
    return(y)
  }

  meanest=apply(est,1,mean)
  biasest=meanest-parameter
  varest=apply(est,1,var)
  mseest=biasest^2+varest
  y=list(meanest=meanest,biasest=biasest,varest=varest,mseest=mseest)
  return(y)
}

# The weighted LS regression file
Weightedls.reg <- function(x,y,W) {
  alpha=0.05

  if(is.data.frame(W)) W=as.vector(t(W))
  if(is.vector(W)) W=diag(W,nrow=length(W))

  m=sum(diag(W))

  if (is.vector(x)){
    n=length(x)
    p=1
    #x<-matrix(x)
    x=t(x)
    x=t(x)

  } else {
    n=dim(x)[1]
    p=dim(x)[2]

    if(p==1) {

      # x<-matrix(x)
      x=x[,1]
      x=t(x)
      x=t(x)

    } else {

      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }
  if (is.vector(y)){

    y=t(y)
    y=t(y)

  } else {
    cy=dim(y)[2]
    if(cy==1) {

      y=y[,1]
      y=t(y)
      y=t(y)

    } else {

      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }

  x=as.matrix(x)
  y=as.matrix(y)

  X <- cbind(matrix(1,n),x)

  Xw <- W^0.5%*%X
  yw <- W^0.5%*%y
  XpX <- t(X)%*%W%*%X
  Xpy <- t(X)%*%W%*%y

  invXpX <- solve(XpX)
  beta <- invXpX%*%Xpy

  yhat <- X%*%beta
  yhatw <- Xw%*%beta

  e <- y-yhat
  ew <- W^(1/2)%*%e
  SSE <- t(ew)%*%ew
  SSE <- as.numeric(SSE)
  MSE <- SSE/(n-(p+1))

  s <- 0
  for (i in 1:n) {
    s <- s+W[i,i]*y[i]
  }

  ymeanw <- s/m

  SST <- t(yw)%*%yw-m*ymeanw^2
  SST <- as.numeric(SST)
  MST <- SST/(n-1)
  SSR <- t(yhatw)%*%yhatw-m*ymeanw^2
  SSR <- as.numeric(SSR)
  R2 <- SSR/SST
  MSR <- SSR/p
  F <- MSR/MSE
  R2adj <- 1-MSE/MST

  sig <- 1-pf(F,p,n-(p+1))
  varbeta <- invXpX*MSE
  stdbeta <- sqrt(diag(varbeta))

  confint <- rbind(t(beta)-qt(1-alpha/2,n-(p+1))*stdbeta,t(beta)+qt(1-alpha/2,n-(p+1))*stdbeta)
  anovatable <- data.frame("s.v."=c("Regression","Error","Total"),
                           "S.S."=c(SSR,SSE,SST),
                           "d.f."=c(p,n-(p+1),n-1),
                           "M.S."=c(MSR,MSE,MST),
                           "F"=c(F,NA,NA),
                           "sig."=c(sig,NA,NA))

  z <- list(beta=beta,e=e,ew=ew,yhat=yhat,yhatw=yhatw,MSE=MSE,F=F,sig=sig,varbeta=varbeta,stdbeta=stdbeta,
            R2=R2,R2adj=R2adj,anovatable=anovatable,confint=confint)
  return(z)
}
# End of weighted LS regression file

# Ridge regression
Ridge.reg <- function(x,y,c) {
  if (is.vector(x)){
    n=length(x)
    p=1
    x=t(x)
    x=t(x)

  } else {
    n=dim(x)[1]
    p=dim(x)[2]

    if(p==1) {

      x=x[,1]
      x=t(x)
      x=t(x)

    } else {

      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }

  if (is.vector(y)){
    y=t(y)
    y=t(y)
  } else {
    cy=dim(y)[2]
    if(cy==1) {

      y=y[,1]
      y=t(y)
      y=t(y)

    } else {

      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }
  x=as.matrix(x)
  y=as.matrix(y)

  yr <- scale(y)/sqrt(n-1)
  xr <- scale(x)/sqrt(n-1)
  X <- xr
  XpX <- t(X)%*%X
  XpXplusc <- XpX+c*diag(p)

  Xpy <- t(X)%*%yr
  p1 <- qr(XpX)$rank-1

  invXpXplusc <- solve(XpXplusc)
  beta <- invXpXplusc%*%Xpy

  tsdx <- apply(x,2,sd)
  betaor <- sd(y)*beta/as.matrix(tsdx)
  beta0 <- c(mean(y),apply(x,2,mean))%*%rbind(1,-betaor)
  betaor <- rbind(beta0,betaor)
  yhator <- cbind(matrix(1,n),x)%*%betaor

  yhat <- X%*%beta
  e <- yr-yhat
  SSE <- t(e)%*%e
  SSE=as.numeric(SSE)
  MSE <- SSE/(n-(p+1))

  s1 <- 0

  for(i in 1:n) {
    s1 <- s1+yr[i]
  }

  ymeanr <- s1/n

  SST <- t(yr)%*%yr-n*ymeanr^2
  SST <- as.numeric(SST)
  MST <- SST/(n-1)

  SSR <- 1-SSE
  SSR <- as.numeric(SSR)
  MSR <- SSR/p
  R2 <- SSR/SST
  R2adj <- 1-MSE/MST

  F <- MSR/MSE
  sig <- 1-pf(F,p,n-(p+1))

  varbeta=invXpXplusc*MSE
  VIF=invXpXplusc%*%XpX%*%invXpXplusc

  aF <- c(F,NA,NA)
  asig <- c(sig,NA,NA)
  s.v <- c("Regression","Error","Total")
  S.S <- c(SSR,SSE,SST)
  d.f <- c(p,n-(p+1),n-1)
  M.S <- c(MSR,MSE,MST)
  anovatable <- data.frame(s.v,S.S,d.f,M.S,aF,asig)

  z <- list(beta=beta,e=e,yhat=yhat,MSE=MSE,F=F,sig=sig,varbeta=varbeta,R2=R2,R2adj=R2adj,anovatable=anovatable,VIF=VIF,betaor=betaor,yhator=yhator)

  return(z)
}
# End of Ridge regression file


# Ridge regression with a selected constant k
ridgereg.k <- function(x,y,a,b) {

  if ((a>b) | (a<0) | (b>1)){
    print('Wrong input, please try again!')
    return(ridgereg.k)

  } else {

    for (j in 1:4) {
      k <- seq(a,b,length=11)
      v <- matrix(0,11,dim(x)[2])

      for (i in 1:11) {

        VIF <- Ridge.reg(x,y,k[i])$VIF
        diagVIF <- as.matrix(diag(VIF))
        tdiagVIF <- t(diagVIF)

        if (min(tdiagVIF)>1) {

          v[i,] <- tdiagVIF[1,]

        } else {

          break
        }
      }
      b <- k[i]
      a<-k[i-1]
      k <- as.matrix(k)
    }

    kk <- k[i-1]
  }

  ridgereg <- Ridge.reg(x,y,k[i-1])
  ridgereg$k=k[i-1]

  z <- list(k=k)
  return(ridgereg)
}
# End of Ridge regression with a selected constant k

#weighted ridge regression file
Weightedridge.reg <- function(x,y,W) {
  alpha=0.05
  a=0
  b=1
  if(is.data.frame(W)) W=as.vector(t(W))
  if(is.vector(W)) W=diag(W,nrow=length(W))

  m=sum(diag(W))

  if (is.vector(x)){
    n=length(x)
    p=1
    x<-matrix(x)

  } else {
    n=dim(x)[1]
    p=dim(x)[2]

    if(p==1) {

      x<-matrix(x)

    } else {

      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }
  if (is.vector(y)){

    y<-matrix(y)

  } else {
    cy=dim(y)[2]
    if(cy==1) {

      y<-matrix(y)

    } else {

      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }

  ridgeregc= ridgereg.k(x,y,a,b)
  cc=ridgeregc$k

  xw <- W^0.5%*%x
  yw <- W^0.5%*%y
  yr <- scale(y)/sqrt(n-1)
  xr <- scale(x)/sqrt(n-1)
  yrw=W^0.5%*%yr
  xrw=W^0.5%*%xr
  X <- xr
  XpX <- t(X)%*%W%*%xr
  XpXplusc <- XpX+cc*diag(p)

  Xpy <- t(X)%*%W%*%yr
  invXpXplusc <- solve(XpXplusc)
  beta <- invXpXplusc%*%Xpy

  tsdx <- apply(xw,2,sd)
  betaor <- sd(yw)*beta/as.matrix(tsdx)
  beta0 <- c(mean(yw),apply(xw,2,mean))%*%rbind(1,-betaor)
  betaor <- rbind(beta0,betaor)
  yhator <- cbind(matrix(1,n),x)%*%betaor

  yhat <- xr%*%beta
  yhatw <- xrw%*%beta

  e <- yr-yhat
  ew <- W^(1/2)%*%e
  SSE <- t(ew)%*%ew
  SSE <- as.numeric(SSE)
  MSE <- SSE/(n-(p+1))

  s <- 0
  for (i in 1:n) {
    s <- s+W[i,i]*yr[i]
  }
  ymeanw <- s/m

  SST <- t(yrw)%*%yrw-m*ymeanw^2
  SST <- as.numeric(SST)
  MST <- SST/(n-1)
  SSR <- SST-SSE
  SSR <- as.numeric(SSR)
  R2 <- SSR/SST
  MSR <- SSR/p
  F <- MSR/MSE
  R2adj <- 1-MSE/MST

  sig <- 1-pf(F,p,n-(p+1))

  varbeta <- invXpXplusc*MSE
  stdbeta <- sqrt(diag(varbeta))
  #print(cc)
  confint <- rbind(t(beta)-qt(1-alpha/2,n-(p+1))*stdbeta,t(beta)+qt(1-alpha/2,n-(p+1))*stdbeta)
  anovatable <- data.frame("s.v."=c("Regression","Error","Total"),
                           "S.S."=c(SSR,SSE,SST),
                           "d.f."=c(p,n-(p+1),n-1),
                           "M.S."=c(MSR,MSE,MST),
                           "F"=c(F,NA,NA),
                           "sig."=c(sig,NA,NA))

  z <- list(cc=cc,beta=beta,betaor=betaor,e=e,ew=ew,yhat=yhat,yhatw=yhatw,yhator=yhator,MSE=MSE,F=F,sig=sig,varbeta=varbeta,stdbeta=stdbeta,
            R2=R2,R2adj=R2adj,anovatable=anovatable,confint=confint)

  return(z)
}
# End of weighted ridge regression file


###-----------------------------------------------------------------------
# Reg S-type estimators program and weighted regression commands end

##the functions needed to run regstype start here
f1=function(u,c) {
  
  2*(((u^2)/2-(u^4)/(2*(c^2))+(u^6)/(6*(c^4)))*((1/sqrt(2*pi)*exp((-u^2)/2))))
  
}

f2=function(u,c) {
  2*((c^2/6)*((1/sqrt(2*pi)*exp((-u^2)/2))))
}
##the functions needed to run regstype end


# Reg S-type estimators program start here
regstype=function(y,x,c) {

  maxit=100
  eps=0.00001
  
  K=integrate(f1,c=c,0,c)$value+integrate(f2,c=c,c,Inf)$value
  
  maxit=100
  eps=0.00001
  
  if (is.vector(x)){
    n=length(x)
    p=1
    x=t(x)
    x=t(x)
    
  } else {
    n=dim(x)[1]
    p=dim(x)[2]
    
    if(p==1) {
      
      x=x[,1]
      x=t(x)
      x=t(x)
      
    } else {
      
      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }
  
  if (is.vector(y)){
    y=t(y)
    y=t(y)
  } else {
    cy=dim(y)[2]
    if(cy==1) {
      
      y=y[,1]
      y=t(y)
      y=t(y)
      
    } else {
      
      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }
  
  s=p+1
  
  regls=lm(y~x)
  
  betals=regls$coefficients
  betals=t(betals)
  betals=t(betals)
  
  els=regls$residuals
  els=t(els)
  els=t(els)
  
  betas=array(NA,dim=c(s,maxit))
  es=array(NA,dim=c(n,maxit))
  us=array(NA,dim=c(n,maxit))
  sigmas=array(NA,dim=c(maxit))
  conds=array(NA,dim=c(maxit))
  
  
  for(i in 1:s){
    betas[i,1]=betals[i]
  }
  
  for(i in 1:n){
    es[i,1]=els[i]
  }
  
  sigmas[1]=mad(es[,1])
  
  for(i in 1:n){
    us[i,1]=es[i,1]/sigmas[1]
  }

  W1s=(1-((us/c)^2))^2
  
  Ws=(abs(us)<=c)*W1s
  
  regtemp=regweighteds(y,x,as.vector(Ws[,1]))
  
  for(i in 1:s){
    betas[i,2]=regtemp$beta[i]
  }
  
  for(i in 1:n){
    es[i,2]=regtemp$e[i]
  }
  
  
  fark=betas[,2]-betas[,1]
  
  conds[1]=norm(t(fark),"2")/norm(t(betas[,2]),"2")
  
  ites=2
  
  while ((conds[ites-1]>=eps)&ites<100) {
    
    sigmas[ites]=sqrt((1/(n*K))*sum(regtemp$ew^2))
    
    for(i in 1:n){
      us[i,ites]=es[i,ites]/sigmas[ites]
      Ws[i,ites]=(((abs(us[i,ites])<=c)*(((us[i,ites]^2)/2)-((us[i,ites]^4)/(2*(c^2)))+((us[i,ites]^6)/(6*(c^4)))))/(us[i,ites]^2))+
        ((abs(us[i,ites])>c)*(((c^2)/6)/(us[i,ites]^2)))
    }
    
    regtemp=regweighteds(y,x,as.vector(Ws[,ites]))
    
    ites=ites+1
    
    for(i in 1:s){
      betas[i,ites]=regtemp$beta[i]
    }
    
    for(i in 1:n){
      es[i,ites]=regtemp$e[i]
    }
    
    fark=betas[,ites]-betas[,ites-1]
    
    conds[ites-1]=norm(t(fark),"2")/norm(t(betas[,ites]),"2")
  }
  
  beta=betas[,ites]
  
  sigma=sigmas[ites-1]
  W=as.vector(Ws[,ites-1])
  
  e=regtemp$e
  yhat=regtemp$yhat
  MSE=regtemp$MSE
  F=regtemp$F
  sig=regtemp$sig
  varbeta=regtemp$varbeta
  stdbeta=regtemp$stdbeta
  R2=regtemp$R2
  R2adj=regtemp$R2adj
  anovatable=regtemp$anovatable
  confint=regtemp$confint
  
  z=list(beta=beta,betas=betas,e=e,es=es,yhat=yhat,MSE=MSE,F=F,sig=sig,varbeta=varbeta,
         stdbeta=stdbeta,R2=R2,R2adj=R2adj,anovatable=anovatable,confint=confint,ites=ites,sigmas=sigmas,
         sigma=sigma,W=W,Ws=Ws,conds=conds)
  
  return(z)
}
# Reg S-type estimators program end


# Weighted regression commands start here
regweighteds=function(y,x,W) {
  
  alpha=0.05
  #alpha=(0.05)/4
  
  if(is.data.frame(W)) W=as.vector(t(W))
  if(is.vector(W)) W=diag(W,nrow=length(W))
  
  m=sum(diag(W))
  
  if (is.vector(x)){
    n=length(x)
    p=1
    x=t(x)
    x=t(x)
    
  } else {
    n=dim(x)[1]
    p=dim(x)[2]
    
    if(p==1) {
      
      x=x[,1]
      x=t(x)
      x=t(x)
      
    } else {
      
      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }
  
  if (is.vector(y)){
    y=t(y)
    y=t(y)
  } else {
    cy=dim(y)[2]
    if(cy==1) {
      
      y=y[,1]
      y=t(y)
      y=t(y)
      
    } else {
      
      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }
  
  
  #  print(x)
  #  print(y)
  
  #print(n)
  #print(p)
  
  #e=sprintf("The current model: Yhat=%f+%f(X-%f) \n \n",
  #          bet[1],bet[2],xM)
  #cat(e)
  
  
  X=cbind(matrix(1,n,1),x)
  
  Xw=W^0.5%*%X
  yw=W^0.5%*%y
  XpX=t(X)%*%W%*%X
  Xpy=t(X)%*%W%*%y
  
  #print(Xw)
  #print(yw)
  #print(XpX)
  #print(Xpy)
  
  invXpX=solve(XpX)
  beta=invXpX%*%Xpy
  
  #print(beta)
  
  yhat=X%*%beta
  yhatw=Xw%*%beta
  
  e=y-yhat
  ew=W^(1/2)%*%e
  
  
  SSE=t(ew)%*%ew
  SSE=as.numeric(SSE)
  MSE=SSE/(n-(p+1))
  
  
  s1=0
  for(i in 1:n) {
    s1=s1+W[i,i]*y[i]
  }
  ymeanw=s1/m
  SST=t(yw)%*%yw-m*ymeanw^2
  SST=as.numeric(SST)
  
  MST=SST/(n-1)
  SSR=t(yhatw)%*%yhatw-m*ymeanw^2
  SSR=as.numeric(SSR)
  
  R2=SSR/SST
  MSR=SSR/p
  F=MSR/MSE
  R2adj=1-MSE/MST
  
  sig=1-pf(F,p,n-(p+1))
  
  varbeta=invXpX*MSE
  stdbeta=sqrt(diag(varbeta))
  confint=rbind(t(beta)-qt(1-alpha/2,n-(p+1))*stdbeta,t(beta)+qt(1-alpha/2,n-(p+1))*stdbeta)
  anovatable=data.frame("s.v."=c("Regression","Error","Total"),
                        "S.S."=c(SSR,SSE,SST),
                        "d.f."=c(p,n-(p+1),n-1),
                        "M.S."=c(MSR,MSE,MST),
                        "F"=c(F,NA,NA),
                        "sig."=c(sig,NA,NA))
  
  z=list(beta=beta,e=e,ew=ew,yhat=yhat,yhatw=yhatw,MSE=MSE,F=F,sig=sig,varbeta=varbeta,
         stdbeta=stdbeta,R2=R2,R2adj=R2adj,anovatable=anovatable,confint=confint)
  
  return(z)
}
# Weighted regression commands end

# Reg S-type estimators program and weighted regression commands end
###-----------------------------------------------------------------------


####robustweightedRidge (the main function to obtain real life data ridge regression results) starts here
robustweightedRidge <- function(x,y) {   
  a=0
  b=1
  
  v <- 32
  k <- 30
  if (is.vector(x)){
    n=length(x)
    p=1
    x=t(x)
    x=t(x)
    
  } else {
    n=dim(x)[1]
    p=dim(x)[2]
    
    if(p==1) {
      x=x[,1]
      x=t(x)
      x=t(x)
    } else {
      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }
  
  if (is.vector(y)){
    y=t(y)
    y=t(y)
    
  } else {
    cy=dim(y)[2]
    if(cy==1) {
      y=y[,1]
      y=t(y)
      y=t(y)
    } else {
      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }
  x=as.matrix(x)
  y=as.matrix(y)
  
  regsr=rlm(y~x,method="MM",maxit=0)
  regsrW=regsr$w
  regsrtemp <- Weightedls.reg(x,y,regsrW)
  e <- regsrtemp$e
  
  stype=regstype(y,x,1.548)
  wstype=stype$W
  weightor=wstype
  weightor=as.vector(weightor)
  WeightedRidgeStype = Weightedridge.reg(x,y,weightor)
  weightorStype=as.matrix(t(t(weightor)))
  e=WeightedRidgeStype$e
  
  regls=Ridge.reg(x,y,0)
                       varbetaols=regls$varbeta
                       stdbetaols <- sqrt(diag(varbetaols))
  
  reglsr=ridgereg.k(x,y,a,b)
  varbetals=reglsr$varbeta
  stdbetals <- sqrt(diag(varbetals))
  
  regtukeymr=rlm(y~x,method="M",psi=psi.bisquare)
  regtukeymrW=regtukeymr$w
  WeightedRidgeTukeym = Weightedridge.reg(x,y,regtukeymrW)
  weightorTukeym=as.matrix(t(t(regtukeymrW)))
  
  reghubermr=rlm(y~x,method="M",psi=psi.huber)
  reghubermrW=reghubermr$w
  WeightedRidgeHuberm = Weightedridge.reg(x,y,reghubermrW)
  weightorHuberm=as.matrix(t(t(reghubermrW)))
  
  regsr=rlm(y~x,method="MM",maxit=0)
  regsrW=regsr$w
  WeightedRidgeS = Weightedridge.reg(x,y,regsrW)
  weightorS=as.matrix(t(t(regsrW)))
  
  regmmr=rlm(y~x,method="MM")
  regmmrW=regmmr$w
  WeightedRidgeMM = Weightedridge.reg(x,y,regmmrW)
  weightorMM=as.matrix(t(t(regmmrW)))
  
  ccStype=WeightedRidgeStype$cc
  ccTukeym=WeightedRidgeTukeym$cc
  ccHuberm=WeightedRidgeHuberm$cc
  ccS=WeightedRidgeS$cc
  ccMM=WeightedRidgeMM$cc
  
                betaols=regls$beta
  betals=reglsr$beta
  betaStype=WeightedRidgeStype$beta
  betaTukeym=WeightedRidgeTukeym$beta
  betaHuberm=WeightedRidgeHuberm$beta
  betaS=WeightedRidgeS$beta
  betaMM=WeightedRidgeMM$beta
  
                betaorols=regls$betaor
  betaorls=reglsr$betaor
  betaorStype=WeightedRidgeStype$betaor
  betaorTukeym=WeightedRidgeTukeym$betaor
  betaorHuberm=WeightedRidgeHuberm$betaor
  betaorS=WeightedRidgeS$betaor
  betaorMM=WeightedRidgeMM$betaor
  
                eols=regls$e
  els=reglsr$e
  eStype=WeightedRidgeStype$e
  eTukeym=WeightedRidgeTukeym$e
  eHuberm=WeightedRidgeHuberm$e
  eS=WeightedRidgeS$e
  eMM=WeightedRidgeMM$e
  
                ewols=regls$ew
  ewls=reglsr$ew
  ewStype=WeightedRidgeStype$ew
  ewTukeym=WeightedRidgeTukeym$ew
  ewHuberm=WeightedRidgeHuberm$ew
  ewS=WeightedRidgeS$ew
  ewMM=WeightedRidgeMM$ew
  
                yhatols=regls$yhat
  yhatls=reglsr$yhat
  yhatStype=WeightedRidgeStype$yhat
  yhatTukeym=WeightedRidgeTukeym$yhat
  yhatHuberm=WeightedRidgeHuberm$yhat
  yhatS=WeightedRidgeS$yhat
  yhatMM=WeightedRidgeMM$yhat
  
                yhatwls=regls$yhatw
  yhatwls=reglsr$yhatw
  yhatwStype=WeightedRidgeStype$yhatw
  yhatwTukeym=WeightedRidgeTukeym$yhatw
  yhatwHuberm=WeightedRidgeHuberm$yhatw
  yhatwS=WeightedRidgeS$yhatw
  yhatwMM=WeightedRidgeMM$yhatw
  
                MSEols=regls$MSE
  MSEls=reglsr$MSE
  MSEStype=WeightedRidgeStype$MSE
  MSETukeym=WeightedRidgeTukeym$MSE
  MSEHuberm=WeightedRidgeHuberm$MSE
  MSES=WeightedRidgeS$MSE
  MSEMM=WeightedRidgeMM$MSE
  
                Fols=regls$F
  Fls=reglsr$F
  FStype=WeightedRidgeStype$F
  FTukeym=WeightedRidgeTukeym$F
  FHuberm=WeightedRidgeHuberm$F
  FS=WeightedRidgeS$F
  FMM=WeightedRidgeMM$F
  
                sigols=regls$sig
  sigls=reglsr$sig
  sigStype=WeightedRidgeStype$sig
  sigTukeym=WeightedRidgeTukeym$sig
  sigHuberm=WeightedRidgeHuberm$sig
  sigS=WeightedRidgeS$sig
  sigMM=WeightedRidgeMM$sig
  
                varbetaols=regls$varbeta
  varbetals=reglsr$varbeta
  varbetaStype=WeightedRidgeStype$varbeta
  varbetaTukeym=WeightedRidgeTukeym$varbeta
  varbetaHuberm=WeightedRidgeHuberm$varbeta
  varbetaS=WeightedRidgeS$varbeta
  varbetaMM=WeightedRidgeMM$varbeta
  
                stdbetaols=stdbetaols
  stdbetals=stdbetals
  stdbetaStype=WeightedRidgeStype$stdbeta
  stdbetaTukeym=WeightedRidgeTukeym$stdbeta
  stdbetaHuberm=WeightedRidgeHuberm$stdbeta
  stdbetaS=WeightedRidgeS$stdbeta
  stdbetaMM=WeightedRidgeMM$stdbeta
  
                R2ols=regls$R2
  R2ls=reglsr$R2
  R2Stype=WeightedRidgeStype$R2
  R2Tukeym=WeightedRidgeTukeym$R2
  R2Huberm=WeightedRidgeHuberm$R2
  R2S=WeightedRidgeS$R2
  R2MM=WeightedRidgeMM$R2
  
                R2adjols=regls$R2adj
  R2adjls=reglsr$R2adj
  R2adjStype=WeightedRidgeStype$R2adj
  R2adjTukeym=WeightedRidgeTukeym$R2adj
  R2adjHuberm=WeightedRidgeHuberm$R2adj
  R2adjS=WeightedRidgeS$R2adj
  R2adjMM=WeightedRidgeMM$R2adj
  
               anovatableols=regls$anovatable
  anovatablels=reglsr$anovatable
  anovatableStype=WeightedRidgeStype$anovatable
  anovatableTukeym=WeightedRidgeTukeym$anovatable
  anovatableHuberm=WeightedRidgeHuberm$anovatable
  anovatableS=WeightedRidgeS$anovatable
  anovatableMM=WeightedRidgeMM$anovatable
  
                confintols=regls$confint
  confintls=reglsr$confint
  confintStype=WeightedRidgeStype$confint
  confintTukeym=WeightedRidgeTukeym$confint
  confintHuberm=WeightedRidgeHuberm$confint
  confintS=WeightedRidgeS$confint
  confintMM=WeightedRidgeMM$confint
  
  s=p+1
  
               estols=regls$betaor
               estols[s+1]=regls$MSE
               estols[s+2]=regls$R2
               estols[s+3]=regls$R2adj
  
  
  estlsr=reglsr$betaor
  estlsr[s+1]=reglsr$MSE
  estlsr[s+2]=reglsr$R2
  estlsr[s+3]=reglsr$R2adj
  
  estStyper=WeightedRidgeStype$betaor
  
  kk=30
  pp=16.5
  
  eStyper=e
  sigma2Styper=((2*pp/kk)*sum(weightorStype*(eStyper^2)))/sum(weightorStype)
  estStyper[s+1]=WeightedRidgeStype$MSE
  estStyper[s+2]=WeightedRidgeStype$R2
  estStyper[s+3]=WeightedRidgeStype$R2adj
  
  esttukeymRidge=WeightedRidgeTukeym$betaor
  esttukeymRidge[s+1]=WeightedRidgeTukeym$MSE
  esttukeymRidge[s+2]=WeightedRidgeTukeym$R2
  esttukeymRidge[s+3]=WeightedRidgeTukeym$R2adj
  
  esthubermRidge=WeightedRidgeHuberm$betaor
  esthubermRidge[s+1]=WeightedRidgeHuberm$MSE
  esthubermRidge[s+2]=WeightedRidgeHuberm$R2
  esthubermRidge[s+3]=WeightedRidgeHuberm$R2adj
  
  estsRidge=WeightedRidgeS$betaor
  estsRidge[s+1]=WeightedRidgeS$MSE
  estsRidge[s+2]=WeightedRidgeS$R2
  estsRidge[s+3]=WeightedRidgeS$R2adj
  
  estmmRidge=WeightedRidgeMM$betaor
  estmmRidge[s+1]=WeightedRidgeMM$MSE
  estmmRidge[s+2]=WeightedRidgeMM$R2
  estmmRidge[s+3]=WeightedRidgeMM$R2adj
  
                  stdols=stdbetaols
  stdlsr=stdbetals
  stdStyper=WeightedRidgeStype$stdbeta
  stdtukeymRidge=WeightedRidgeTukeym$stdbeta
  stdhubermRidge=WeightedRidgeHuberm$stdbeta
  stdsRidge=WeightedRidgeS$stdbeta
  stdmmRidge=WeightedRidgeMM$stdbeta
 
  esttable=cbind(estols,estlsr,estStyper,esthubermRidge, esttukeymRidge, estsRidge,estmmRidge)
  esttabledataframe=data.frame(esttable)
  
  colnames(esttabledataframe)=c("LS","LSridge","StypeR","Huber MRidge","Tukey MRidge","SRidge","MMRidge")
  
  # add row names
  row_names <- c("Intercept", paste0("X.", 1:p), "MSE", "R2", "adj R2")
  rownames(esttabledataframe) <- row_names
  
  stdtable=cbind(stdols,stdlsr,stdStyper,stdhubermRidge,stdtukeymRidge, stdsRidge,stdmmRidge)
  stdtabledataframe=data.frame(stdtable)
  
  colnames(stdtabledataframe)=c("LS","LSridge","StypeR","Huber MRidge","Tukey MRidge","SRidge","MMRidge")
  # add row names in stdtable
  rownames(stdtabledataframe) <- paste0("X", 1:p)
  
  z=list(esttabledataframe=esttabledataframe,
         stdtabledataframe=stdtabledataframe,
         betaols=betaols,betaorols=betaorols,eols=eols,ewols=ewols,yhatols=yhatols,
        MSEols=MSEols,Fols=Fols,sigols=sigols,varbetaols=varbetaols,
         stdbetaols=stdbetaols,R2ols=R2ols,R2adjols=R2adjols,anovatableols=anovatableols,
         confintols=confintols,
         betals=betals,betaorls=betaorls,els=els,ewls=ewls,yhatls=yhatls,
         yhatwls=yhatwls,MSEls=MSEls,Fls=Fls,sigls=sigls,varbetals=varbetals,
         stdbetals=stdbetals,R2ls=R2ls,R2adjls=R2adjls,anovatablels=anovatablels,
         confintls=confintls,
         ccStype=ccStype,betaStype=betaStype,betaorStype=betaorStype,eStype=eStype,ewStype=ewStype,yhatStype=yhatStype,
         yhatwStype=yhatwStype,MSEStype=MSEStype,FStype=FStype,sigStype=sigStype,varbetaStype=varbetaStype,
         stdbetaStype=stdbetaStype,R2Stype=R2Stype,R2adjStype=R2adjStype,anovatableStype=anovatableStype,
         confintStype=confintStype,weightorStype=weightorStype,
         ccTukeym=ccTukeym,betaTukeym=betaTukeym,betaorTukeym=betaorTukeym,eTukeym=eTukeym,ewTukeym=ewTukeym,yhatTukeym=yhatTukeym,
         yhatwTukeym=yhatwTukeym,MSETukeym=MSETukeym,FTukeym=FTukeym,sigTukeym=sigTukeym,varbetaTukeym=varbetaTukeym,
         stdbetaTukeym=stdbetaTukeym,R2Tukeym=R2Tukeym,R2adjTukeym=R2adjTukeym,anovatableTukeym=anovatableTukeym,
         confintTukeym=confintTukeym,weightorTukeym=weightorTukeym,
         ccHuberm=ccHuberm,betaHuberm=betaHuberm,betaorHuberm=betaorHuberm,eHuberm=eHuberm,ewHuberm=ewHuberm,yhatHuberm=yhatHuberm,
         yhatwHuberm=yhatwHuberm,MSEHuberm=MSEHuberm,FHuberm=FHuberm,sigHuberm=sigHuberm,varbetaHuberm=varbetaHuberm,
         stdbetaHuberm=stdbetaHuberm,R2Huberm=R2Huberm,R2adjHuberm=R2adjHuberm,anovatableHuberm=anovatableHuberm,
         confintHuberm=confintHuberm,weightorHuberm=weightorHuberm,
         ccS=ccS,betaS=betaS,betaorS=betaorS,Se=eS,ewS=ewS,yhatS=yhatS,yhatwS=yhatwS,MSES=MSES,FS=FS,sigS=sigS,varbetaS=varbetaS,
         stdbetaS=stdbetaS,R2S=R2S,R2adjS=R2adjS,anovatableS=anovatableS,confintS=confintS,weightorS=weightorS,
         ccMM=ccMM,betaMM=betaMM,betaorMM=betaorMM,eMM=eMM,ewMM=ewMM,yhatMM=yhatMM,yhatwMM=yhatwMM,MSEMM=MSEMM,FMM=FMM,sigMM=sigMM,
         varbetaMM=varbetaMM,stdbetaMM=stdbetaMM,R2MM=R2MM,R2adjMM=R2adjMM,anovatableMM=anovatableMM,confintMM=confintMM,
         weightorMM=weightorMM)
  
  
  return(z)
  
}
####robustweightedRidge end


# The LS regression file
regmy=function(y,x) {

  alpha=0.05

  if (is.vector(x)){
    n=length(x)
    p=1
    x=t(x)
    x=t(x)

  } else {
    n=dim(x)[1]
    p=dim(x)[2]

    if(p==1) {

      x=x[,1]
      x=t(x)
      x=t(x)

    } else {

      colind=2
      xx=cbind(x[,1],x[,2])
      while(colind<p) {
        colind=colind+1
        xx=cbind(xx,x[,colind])
      }
      x=xx }
  }

  if (is.vector(y)){
    y=t(y)
    y=t(y)
  } else {
    cy=dim(y)[2]
    if(cy==1) {

      y=y[,1]
      y=t(y)
      y=t(y)

    } else {

      colind=2
      yy=cbind(y[,1],y[,2])
      while(colind<cy) {
        colind=colind+1
        yy=cbind(yy,y[,colind])
      }
      y=yy }
  }


  X=cbind(matrix(1,n,1),x)

  XpX=t(X)%*%X
  Xpy=t(X)%*%y

  invXpX=solve(XpX)
  beta=invXpX%*%Xpy

  yhat=X%*%beta

  e=y-yhat

  SSE=t(e)%*%e
  SSE=as.numeric(SSE)
  MSE=SSE/(n-(p+1))

  s1=0
  for(i in 1:n) {
    s1=s1+y[i]
  }

  ymean=s1/n
  SST=t(y)%*%y-n*ymean^2
  SST=as.numeric(SST)

  MST=SST/(n-1)
  SSR=t(yhat)%*%yhat-n*ymean^2
  SSR=as.numeric(SSR)

  R2=SSR/SST
  MSR=SSR/p
  F=MSR/MSE
  R2adj=1-MSE/MST

  sig=1-pf(F,p,n-(p+1))

  varbeta=invXpX*MSE
  stdbeta=sqrt(diag(varbeta))
  confint=rbind(t(beta)-qt(1-alpha/2,n-(p+1))*stdbeta,t(beta)+qt(1-alpha/2,n-(p+1))*stdbeta)
  anovatable=data.frame("s.v."=c("Regression","Error","Total"),
                        "S.S."=c(SSR,SSE,SST),
                        "d.f."=c(p,n-(p+1),n-1),
                        "M.S."=c(MSR,MSE,MST),
                        "F"=c(F,NA,NA),
                        "sig."=c(sig,NA,NA))

  z=list(beta=beta,e=e,yhat=yhat,MSE=MSE,F=F,sig=sig,varbeta=varbeta,
         stdbeta=stdbeta,R2=R2,R2adj=R2adj,anovatable=anovatable,confint=confint)

  return(z)
}
# The LS regression file end


##ordermultiplewithranks start here
ordermultiplewithranks=function(x){
  
  if (is.vector(x)){
    
    xs=sort(x)
    n=length(x)
    rankx=array(NA,dim=c(1,n))
    for(i in 1:n){
      #rankx[xs[i]==x]=i
      rankx[order(x)[i]]=i
    }
    
    z=list(xs=xs,rankx=rankx)
    return(z)
  }
  
  dimx=dim(x)
  n=dimx[1]
  nn=dimx[2]
  xs=array(0,dim=c(n,nn))
  rankx=array(NA,dim=c(n,nn))
  
  for(j in 1:nn) {
    
    for(i in 1:n){
      
      rankx[order(x[,j])[i],j]=i
      xs[,j]=x[order(x[,j]),j]
    }
    
  }
  
  z=list(xs=xs,rankx=rankx)
  return(z)
}
##ordermultiplewithranks end


# Ridge Regression Simulation with Robust Methods starts here
ridgeglr.simu=function(n,p,nn) {
  a=0
  b=1
  s=(p+1)
  
  k=3
  pi=0.1
  
  MSE=1
  stdbeta0=1
  stdbeta1=1
  stdbeta2=1
  stdbeta3=1
  
  kk=30
  pp=16.5
  
  # x=rnorm(n*nn)
  # x=array(x,dim=c(n,p,nn))

  mu=rep(0,times=p)
  v=rep(1,times=p)
  
  #cor=matrix(c(1,0,0,0,1,0,0,0,1),nrow = p)
  #cor=matrix(c(1,0.50,0.50,0.50,1,0.50,0.50,0.50,1),nrow = p)
  #cor=matrix(c(1,0.80,0.80,0.80,1,0.80,0.80,0.80,1),nrow = p)
  #cor=matrix(c(1,0.95,0.95,0.95,1,0.95,0.95,0.95,1),nrow = p)
  #cor=matrix(c(1,0.98,0.98,0.98,1,0.98,0.98,0.98,1),nrow = p)
  cor=matrix(c(1,0.999,0.999,0.999,1,0.999,0.999,0.999,1),nrow = p)
  #cor=matrix(c(1,0.99,0.99,0.99,1,0.99,0.99,0.99,1),nrow = p)
  #cor=matrix(c(1,0.75,0.85,0.75,1,0.98,0.85,0.98,1),nrow = p)
  
  mysigma=(diag(v))^(1/2)%*%cor%*%(diag(v))^(1/2)
  
  # for uniform x
  #x=runif(n*p*nn)
  #x=array(x,dim=c(n,p,nn))
  
  # mixture model
  
  #u=runif(n*nn)
  #u=array(u,dim=c(n,1,nn))
  #e3=k*e
  #e=(u<pi)*e3+(u>=pi)*e
  #e=e/sqrt((1-pi)*(1^2)+pi*(k^2))
  
  beta0=as.vector(0)
  beta1=rep(1,times=p)
  beta=c(beta0,beta1)
  
  y=array(NA,dim=c(n,1,nn))
  x=array(NA,dim=c(n,p,nn))
  e=array(NA,dim=c(n,1,nn))
  
  betals=array(NA,dim=c(s,nn))
  betalsridge=array(NA,dim=c(s,nn))
  betaStyper=array(NA,dim=c(s,nn))
  betahubermridge=array(NA,dim=c(s,nn))
  betatukeymridge=array(NA,dim=c(s,nn))
  betasridge=array(NA,dim=c(s,nn))
  betammridge=array(NA,dim=c(s,nn))
  
  MSEls=array(NA,dim=c(1,nn))
  stdbeta1ls=array(NA,dim=c(1,nn))
  stdbeta2ls=array(NA,dim=c(1,nn))
  stdbeta3ls=array(NA,dim=c(1,nn))
  
  MSElsridge=array(NA,dim=c(1,nn))
  stdbeta1lsridge=array(NA,dim=c(1,nn))
  stdbeta2lsridge=array(NA,dim=c(1,nn))
  stdbeta3lsridge=array(NA,dim=c(1,nn))
  
  MSEStype=array(NA,dim=c(1,nn))
  stdbeta1Stype=array(NA,dim=c(1,nn))
  stdbeta2Stype=array(NA,dim=c(1,nn))
  stdbeta3Stype=array(NA,dim=c(1,nn))
  
  MSEhubermr=array(NA,dim=c(1,nn))
  stdbeta1hubermr=array(NA,dim=c(1,nn))
  stdbeta2hubermr=array(NA,dim=c(1,nn))
  stdbeta3hubermr=array(NA,dim=c(1,nn))
  
  MSEtukeymr=array(NA,dim=c(1,nn))
  stdbeta1tukeymr=array(NA,dim=c(1,nn))
  stdbeta2tukeymr=array(NA,dim=c(1,nn))
  stdbeta3tukeymr=array(NA,dim=c(1,nn))
  
  MSEsr=array(NA,dim=c(1,nn))
  stdbeta1sr=array(NA,dim=c(1,nn))
  stdbeta2sr=array(NA,dim=c(1,nn))
  stdbeta3sr=array(NA,dim=c(1,nn))
  
  MSEmmr=array(NA,dim=c(1,nn))
  stdbeta1mmr=array(NA,dim=c(1,nn))
  stdbeta2mmr=array(NA,dim=c(1,nn))
  stdbeta3mmr=array(NA,dim=c(1,nn))

  for(i in 1:nn) {
    
    x[,,i]=mvrnorm(n=n,mu=mu,Sigma=mysigma)
    e[,,i]=rnorm(n)
    
    #outlier model
    r=floor(0.5+pi*n)
    e[1:r,,i]=k*e[1:r,,i]
    e[,,i]=e[,,i]/sqrt((1-(r/n))*(1^2)+(r/n)*(k^2))
    
    y[,,i]=cbind(rep(c(1),times=n),x[,,i])%*%t(t(beta))+e[,,i]
    
    regls=regmy(y[,,i],x[,,i])
    betals[,i]=regls$beta
    
    MSEls[,i]=regls$MSE
    stdbeta1ls[,i]=regls$stdbeta[2]
    stdbeta2ls[,i]=regls$stdbeta[3]
    stdbeta3ls[,i]=regls$stdbeta[4]
    
    regStyper<-robustweightedRidge(x[,,i],y[,,i])
    
    MSElsridge[,i]=regStyper$esttabledataframe$LSridge[5]
    MSEStype[,i]=regStyper$esttabledataframe$StypeR[5]
    MSEhubermr[,i]=regStyper$esttabledataframe$`Huber MRidge`[5]
    MSEtukeymr[,i]=regStyper$esttabledataframe$`Tukey MRidge`[5]
    MSEsr[,i]=regStyper$esttabledataframe$SRidge[5]
    MSEmmr[,i]=regStyper$esttabledataframe$MMRidge[5]
    
    
    betalsridge[,i]= regStyper$betaorls
    betaStyper[,i]=regStyper$betaorStype
    
    stdbeta1lsridge[,i]=regStyper$stdtabledataframe$LSridge[1]
    stdbeta2lsridge[,i]=regStyper$stdtabledataframe$LSridge[2]
    stdbeta3lsridge[,i]=regStyper$stdtabledataframe$LSridge[3]
    
    reghubermr=rlm(y[,,i]~x[,,i],method="M",psi=psi.huber)
    reghubermrW=reghubermr$w
    WeightedRidgeHuberm = Weightedridge.reg(x[,,i],y[,,i],reghubermrW)
    betahubermridge[,i]=WeightedRidgeHuberm$betaor
    
    regtukeymr=rlm(y[,,i]~x[,,i],method="M",psi=psi.bisquare)
    regtukeymrW=regtukeymr$w
    WeightedRidgeTukeym = Weightedridge.reg(x[,,i],y[,,i],regtukeymrW)
    betatukeymridge[,i]=WeightedRidgeTukeym$betaor
    
    regsr=rlm(y[,,i]~x[,,i],method="MM",maxit=0)
    regsrW=regsr$w
    WeightedRidgeS = Weightedridge.reg(x[,,i],y[,,i],regsrW)
    betasridge[,i]=WeightedRidgeS$betaor
    
    regmmr=rlm(y[,,i]~x[,,i],method="MM")
    regmmrW=regmmr$w
    WeightedRidgeMM = Weightedridge.reg(x[,,i],y[,,i],regmmrW)
    betammridge[,i]=WeightedRidgeMM$ betaor
    
    stdbeta1Stype[,i]=regStyper$stdtabledataframe$StypeR[1]
    stdbeta2Stype[,i]=regStyper$stdtabledataframe$StypeR[2]
    stdbeta3Stype[,i]=regStyper$stdtabledataframe$StypeR[3]
    
    stdbeta1hubermr[,i]=regStyper$stdtabledataframe$`Huber MRidge`[1]
    stdbeta2hubermr[,i]=regStyper$stdtabledataframe$`Huber MRidge`[2]
    stdbeta3hubermr[,i]=regStyper$stdtabledataframe$`Huber MRidge`[3]
    
    stdbeta1tukeymr[,i]=regStyper$stdtabledataframe$`Tukey MRidge`[1]
    stdbeta2tukeymr[,i]=regStyper$stdtabledataframe$`Tukey MRidge`[2]
    stdbeta3tukeymr[,i]=regStyper$stdtabledataframe$`Tukey MRidge`[3]
    
    stdbeta1sr[,i]=regStyper$stdtabledataframe$SRidge[1]
    stdbeta2sr[,i]=regStyper$stdtabledataframe$SRidge[2]
    stdbeta3sr[,i]=regStyper$stdtabledataframe$SRidge[3]
    
    stdbeta1mmr[,i]=regStyper$stdtabledataframe$MMRidge[1]
    stdbeta2mmr[,i]=regStyper$stdtabledataframe$MMRidge[2]
    stdbeta3mmr[,i]=regStyper$stdtabledataframe$MMRidge[3]
    
  }
  
  MSElsv=MSEls[1,]
  MSElsridgev=MSElsridge[1,]
  MSEStypev=MSEStype[1,]
  MSEhubermrv=MSEhubermr[1,]
  MSEtukeymrv=MSEtukeymr[1,]
  MSEsrv=MSEsr[1,]
  MSEmmrv=MSEmmr[1,]
  
  stdbeta1lsv=stdbeta1ls[1,]
  stdbeta2lsv=stdbeta2ls[1,]
  stdbeta3lsv=stdbeta3ls[1,]
  
  stdbeta1lsridgev=stdbeta1lsridge[1,]
  stdbeta2lsridgev=stdbeta2lsridge[1,]
  stdbeta3lsridgev=stdbeta3lsridge[1,]
  
  stdbeta1Stypev=stdbeta1Stype[1,]
  stdbeta2Stypev=stdbeta2Stype[1,]
  stdbeta3Stypev=stdbeta3Stype[1,]
  
  stdbeta1hubermrv=stdbeta1hubermr[1,]
  stdbeta2hubermrv=stdbeta2hubermr[1,]
  stdbeta3hubermrv=stdbeta3hubermr[1,]
  
  stdbeta1tukeymrv=stdbeta1tukeymr[1,]
  stdbeta2tukeymrv=stdbeta2tukeymr[1,]
  stdbeta3tukeymrv=stdbeta3tukeymr[1,]
  
  stdbeta1srv=stdbeta1sr[1,]
  stdbeta2srv=stdbeta2sr[1,]
  stdbeta3srv=stdbeta3sr[1,]
  
  stdbeta1mmrv=stdbeta1mmr[1,]
  stdbeta2mmrv=stdbeta2mmr[1,]
  stdbeta3mmrv=stdbeta3mmr[1,]
  
  MSElsv=MSEls[1,]
  MSElsridgev=MSElsridge[1,]
  MSEStypev=MSEStype[1,]
  MSEhubermrv=MSEhubermr[1,]
  MSEtukeymrv=MSEtukeymr[1,]
  MSEsrv=MSEsr[1,]
  MSEmmrv=MSEmmr[1,]
  
  MSElist=rbind(MSElsv,MSElsridgev,MSEStypev,MSEhubermrv,MSEtukeymrv,
                MSEsrv,MSEmmrv)
  
  MSEtemp=ordermultiplewithranks(MSElist)
  
  orderedMSE=MSEtemp$xs
  rankMSE=MSEtemp$rankx
  
  
  stdbeta1list=rbind(stdbeta1lsv,stdbeta1lsridgev,stdbeta1Stypev,stdbeta1hubermrv,stdbeta1tukeymrv,
                     stdbeta1srv,stdbeta1mmrv)
  
  stdbeta1temp=ordermultiplewithranks(stdbeta1list)
  
  orderedstdbeta1=stdbeta1temp$xs
  rankstdbeta1=stdbeta1temp$rankx
  
  
  stdbeta2list=rbind(stdbeta2lsv,stdbeta2lsridgev,stdbeta2Stypev,stdbeta2hubermrv,stdbeta2tukeymrv,
                     stdbeta2srv,stdbeta2mmrv)
  
  stdbeta2temp=ordermultiplewithranks(stdbeta2list)
  
  orderedstdbeta2=stdbeta2temp$xs
  rankstdbeta2=stdbeta2temp$rankx
  
  stdbeta3list=rbind(stdbeta3lsv,stdbeta3lsridgev,stdbeta3Stypev,stdbeta3hubermrv,stdbeta3tukeymrv,
                     stdbeta3srv,stdbeta3mmrv)
  
  stdbeta3temp=ordermultiplewithranks(stdbeta3list)
  
  orderedstdbeta3=stdbeta3temp$xs
  rankstdbeta3=stdbeta3temp$rankx
  
  meanrankMSE=apply(rankMSE,1,mean)
  meanrankMSE=t(meanrankMSE)
  meanrankMSE=t(meanrankMSE)
  
  rownames(meanrankMSE)=c("ls","lsridge","Stype","hubermr","tukeymr","sr","mmr")
  meanrankMSEdataframe=data.frame(meanrankMSE)
  
  meanrankstdbeta1=apply(rankstdbeta1,1,mean)
  meanrankstdbeta1=t(meanrankstdbeta1)
  meanrankstdbeta1=t(meanrankstdbeta1)
  
  meanrankstdbeta1dataframe=data.frame(meanrankstdbeta1)
  
  meanrankstdbeta2=apply(rankstdbeta2,1,mean)
  meanrankstdbeta2=t(meanrankstdbeta2)
  meanrankstdbeta2=t(meanrankstdbeta2)
  
  meanrankstdbeta2dataframe=data.frame(meanrankstdbeta2)
  
  meanrankstdbeta3=apply(rankstdbeta3,1,mean)
  meanrankstdbeta3=t(meanrankstdbeta3)
  meanrankstdbeta3=t(meanrankstdbeta3)
  
  meanrankstdbeta3dataframe=data.frame(meanrankstdbeta3)
  
  
  meanrankdataframe=cbind(meanrankMSEdataframe,meanrankstdbeta1dataframe,meanrankstdbeta2dataframe,
                          meanrankstdbeta3dataframe)
  
  simubetals=meanbiasvarmse(betals,beta)
  simubetalsridge=meanbiasvarmse(betalsridge,beta)
  simubetaStype=meanbiasvarmse(betaStyper,beta)
  simubetahubermr=meanbiasvarmse(betahubermridge,beta)
  simubetatukeymr=meanbiasvarmse(betatukeymridge,beta)
  simubetasr=meanbiasvarmse(betasridge,beta)
  simubetammr=meanbiasvarmse(betammridge,beta)
  
  simuMSEls=meanbiasvarmse(MSEls,MSE)
  simuMSElsridge=meanbiasvarmse(MSElsridge,MSE)
  simuMSEStype=meanbiasvarmse(MSEStype,MSE)
  simuMSEhubermr=meanbiasvarmse(MSEhubermr,MSE)
  simuMSEtukeymr=meanbiasvarmse(MSEtukeymr,MSE)
  simuMSEsr=meanbiasvarmse(MSEsr,MSE)
  simuMSEmmr=meanbiasvarmse(MSEmmr,MSE)
  
  simustdbeta1ls=meanbiasvarmse(stdbeta1ls,stdbeta1)
  simustdbeta2ls=meanbiasvarmse(stdbeta2ls,stdbeta2)
  simustdbeta3ls=meanbiasvarmse(stdbeta3ls,stdbeta3)
  
  simustdbeta1lsridge=meanbiasvarmse(stdbeta1lsridge,stdbeta1)
  simustdbeta2lsridge=meanbiasvarmse(stdbeta2lsridge,stdbeta2)
  simustdbeta3lsridge=meanbiasvarmse(stdbeta3lsridge,stdbeta3)
  
  simustdbeta1Stype=meanbiasvarmse(stdbeta1Stype,stdbeta1)
  simustdbeta2Stype=meanbiasvarmse(stdbeta2Stype,stdbeta2)
  simustdbeta3Stype=meanbiasvarmse(stdbeta3Stype,stdbeta3)
  
  simustdbeta1hubermr=meanbiasvarmse(stdbeta1hubermr,stdbeta1)
  simustdbeta2hubermr=meanbiasvarmse(stdbeta2hubermr,stdbeta2)
  simustdbeta3hubermr=meanbiasvarmse(stdbeta3hubermr,stdbeta3)
  
  simustdbeta1tukeymr=meanbiasvarmse(stdbeta1tukeymr,stdbeta1)
  simustdbeta2tukeymr=meanbiasvarmse(stdbeta2tukeymr,stdbeta2)
  simustdbeta3tukeymr=meanbiasvarmse(stdbeta3tukeymr,stdbeta3)
  
  simustdbeta1sr=meanbiasvarmse(stdbeta1sr,stdbeta1)
  simustdbeta2sr=meanbiasvarmse(stdbeta2sr,stdbeta2)
  simustdbeta3sr=meanbiasvarmse(stdbeta3sr,stdbeta3)
  
  simustdbeta1mmr=meanbiasvarmse(stdbeta1mmr,stdbeta1)
  simustdbeta2mmr=meanbiasvarmse(stdbeta2mmr,stdbeta2)
  simustdbeta3mmr=meanbiasvarmse(stdbeta3mmr,stdbeta3)
  
  simumeanls=simubetals$meanest
  simumeanlsridge=simubetalsridge$meanest
  simumeanStype=simubetaStype$meanest
  simumeanhubermr=simubetahubermr$meanest
  simumeantukeymr=simubetatukeymr$meanest
  simumeansr=simubetasr$meanest
  simumeanmmr=simubetammr$meanest
  
  simubiasls=simubetals$biasest
  simubiaslsridge=simubetalsridge$biasest
  simubiasStype=simubetaStype$biasest
  simubiashubermr=simubetahubermr$biasest
  simubiastukeymr=simubetatukeymr$biasest
  simubiassr=simubetasr$biasest
  simubiasmmr=simubetammr$biasest
  
  simuvarls=n*simubetals$varest
  simuvarlsridge=n*simubetalsridge$varest
  simuvarStype=n*simubetaStype$varest
  simuvarhubermr=n*simubetahubermr$varest
  simuvartukeymr=n*simubetatukeymr$varest
  simuvarsr=n*simubetasr$varest
  simuvarmmr=n*simubetammr$varest
  
  simumsels=n*simubetals$mseest
  simumselsridge=n*simubetalsridge$mseest
  simumseStype=n*simubetaStype$mseest
  simumsehubermr=n*simubetahubermr$mseest
  simumsetukeymr=n*simubetatukeymr$mseest
  simumsesr=n*simubetasr$mseest
  simumsemmr=n*simubetammr$mseest
  
  simumeanls[5]=simuMSEls$meanest
  simumeanls[6]=simustdbeta1ls$meanest
  simumeanls[7]=simustdbeta2ls$meanest
  simumeanls[8]=simustdbeta3ls$meanest
  
  simumeanlsridge[5]=simuMSElsridge$meanest
  simumeanlsridge[6]=simustdbeta1lsridge$meanest
  simumeanlsridge[7]=simustdbeta2lsridge$meanest
  simumeanlsridge[8]=simustdbeta3lsridge$meanest
  
  simumeanStype[5]=simuMSEStype$meanest
  simumeanStype[6]=simustdbeta1Stype$meanest
  simumeanStype[7]=simustdbeta2Stype$meanest
  simumeanStype[8]=simustdbeta3Stype$meanest
  
  simumeanhubermr[5]=simuMSEhubermr$meanest
  simumeanhubermr[6]=simustdbeta1hubermr$meanest
  simumeanhubermr[7]=simustdbeta2hubermr$meanest
  simumeanhubermr[8]=simustdbeta3hubermr$meanest
  
  simumeantukeymr[5]=simuMSEtukeymr$meanest
  simumeantukeymr[6]=simustdbeta1tukeymr$meanest
  simumeantukeymr[7]=simustdbeta2tukeymr$meanest
  simumeantukeymr[8]=simustdbeta3tukeymr$meanest
  
  simumeansr[5]=simuMSEsr$meanest
  simumeansr[6]=simustdbeta1sr$meanest
  simumeansr[7]=simustdbeta2sr$meanest
  simumeansr[8]=simustdbeta3sr$meanest
  
  simumeanmmr[5]=simuMSEmmr$meanest
  simumeanmmr[6]=simustdbeta1mmr$meanest
  simumeanmmr[7]=simustdbeta2mmr$meanest
  simumeanmmr[8]=simustdbeta3mmr$meanest
  
  simubiasls[5]=NA
  simubiasls[6]=NA
  simubiasls[7]=NA
  simubiasls[8]=NA
  
  simubiaslsridge[5]=NA
  simubiaslsridge[6]=NA
  simubiaslsridge[7]=NA
  simubiaslsridge[8]=NA
  
  simubiasStype[5]=NA
  simubiasStype[6]=NA
  simubiasStype[7]=NA
  simubiasStype[8]=NA
  
  simubiashubermr[5]=NA
  simubiashubermr[6]=NA
  simubiashubermr[7]=NA
  simubiashubermr[8]=NA
  
  simubiastukeymr[5]=NA
  simubiastukeymr[6]=NA
  simubiastukeymr[7]=NA
  simubiastukeymr[8]=NA
  
  simubiassr[5]=NA
  simubiassr[6]=NA
  simubiassr[7]=NA
  simubiassr[8]=NA
  
  simubiasmmr[5]=NA
  simubiasmmr[6]=NA
  simubiasmmr[7]=NA
  simubiasmmr[8]=NA
  
  simuvarls[5]=n*simuMSEls$varest
  simuvarls[6]=n*simustdbeta1ls$varest
  simuvarls[7]=n*simustdbeta2ls$varest
  simuvarls[8]=n*simustdbeta3ls$varest
  
  simuvarlsridge[5]=n*simuMSElsridge$varest
  simuvarlsridge[6]=n*simustdbeta1lsridge$varest
  simuvarlsridge[7]=n*simustdbeta2lsridge$varest
  simuvarlsridge[8]=n*simustdbeta3lsridge$varest
  
  simuvarStype[5]=n*simuMSEStype$varest
  simuvarStype[6]=n*simustdbeta1Stype$varest
  simuvarStype[7]=n*simustdbeta2Stype$varest
  simuvarStype[8]=n*simustdbeta3Stype$varest
  
  simuvarhubermr[5]=n*simuMSEhubermr$varest
  simuvarhubermr[6]=n*simustdbeta1hubermr$varest
  simuvarhubermr[7]=n*simustdbeta2hubermr$varest
  simuvarhubermr[8]=n*simustdbeta3hubermr$varest
  
  simuvartukeymr[5]=n*simuMSEtukeymr$varest
  simuvartukeymr[6]=n*simustdbeta1tukeymr$varest
  simuvartukeymr[7]=n*simustdbeta2tukeymr$varest
  simuvartukeymr[8]=n*simustdbeta3tukeymr$varest
  
  simuvarsr[5]=n*simuMSEsr$varest
  simuvarsr[6]=n*simustdbeta1sr$varest
  simuvarsr[7]=n*simustdbeta2sr$varest
  simuvarsr[8]=n*simustdbeta3sr$varest
  
  simuvarmmr[5]=n*simuMSEmmr$varest
  simuvarmmr[6]=n*simustdbeta1mmr$varest
  simuvarmmr[7]=n*simustdbeta2mmr$varest
  simuvarmmr[8]=n*simustdbeta3mmr$varest
  
  simumsels[5]=NA
  simumsels[6]=NA
  simumsels[7]=NA
  simumsels[8]=NA
  
  simumselsridge[5]=NA
  simumselsridge[6]=NA
  simumselsridge[7]=NA
  simumselsridge[8]=NA
  
  simumseStype[5]=NA
  simumseStype[6]=NA
  simumseStype[7]=NA
  simumseStype[8]=NA
  
  simumsehubermr[5]=NA
  simumsehubermr[6]=NA
  simumsehubermr[7]=NA
  simumsehubermr[8]=NA
  
  simumsetukeymr[5]=NA
  simumsetukeymr[6]=NA
  simumsetukeymr[7]=NA
  simumsetukeymr[8]=NA
  
  simumsesr[5]=NA
  simumsesr[6]=NA
  simumsesr[7]=NA
  simumsesr[8]=NA
  
  simumsemmr[5]=NA
  simumsemmr[6]=NA
  simumsemmr[7]=NA
  simumsemmr[8]=NA
  
  simuREffls=100*simumselsridge/simumsels
  simuREfflsridge=100*simumselsridge/simumselsridge
  simuREffStype=100*simumselsridge/simumseStype
  simuREfftukeymr=100*simumselsridge/simumsetukeymr
  simuREffhubermr=100*simumselsridge/simumsehubermr
  simuREffsr=100*simumselsridge/simumsesr
  simuREffmmr=100*simumselsridge/simumsemmr
  
  for(i in 5:8){
    
    simuREffls[i]=100*simumeanlsridge[i]/simumeanls[i]
    simuREfflsridge[i]=100*simumeanlsridge[i]/simumeanlsridge[i]
    simuREffStype[i]=100*simumeanlsridge[i]/simumeanStype[i]
    simuREfftukeymr[i]=100*simumeanlsridge[i]/simumeantukeymr[i]
    simuREffhubermr[i]=100*simumeanlsridge[i]/simumeanhubermr[i]
    simuREffsr[i]=100*simumeanlsridge[i]/simumeansr[i]
    simuREffmmr[i]=100*simumeanlsridge[i]/simumeanmmr[i]
    
  }
  
  simutable=rbind(simumeanls,simumeanlsridge, simumeanStype,simumeanhubermr,simumeantukeymr,
                  simumeansr,simumeanmmr)
  
  simutable=rbind(simutable,simubiasls,simubiaslsridge, simubiasStype,simubiashubermr,simubiastukeymr,
                  simubiassr,simubiasmmr)
  
  simutable=rbind(simutable,simuvarls,simuvarlsridge, simuvarStype,simuvarhubermr,simuvartukeymr,
                  simuvarsr,simuvarmmr)
  
  simutable=rbind(simutable,simumsels,simumselsridge,simumseStype,simumsehubermr,simumsetukeymr,
                  simumsesr,simumsemmr)
  
  simutable=rbind(simutable,simuREffls,simuREfflsridge,simuREffStype,simuREffhubermr,simuREfftukeymr,
                  simuREffsr,simuREffmmr)
  
  simutabledataframe=data.frame(simutable)
  colnames(simutabledataframe)=c("b0","b1","b2","b3","MSE","stdbeta1","stdbeta2","stdbeta3")
 # print(simutabledataframe)
#  print(meanrankdataframe)
  z=list(simutabledataframe=simutabledataframe,meanrankdataframe=meanrankdataframe)
  
}
# Ridge Regression Simulation with Robust Methods end


