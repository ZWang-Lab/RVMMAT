#'  NULL model estimation
#'
#' This function estimate the parameters and residuals for the NULL model in RVMMAT
#'
#' @param y.long Long-formatted phenotype vector
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param phe.model String, the phenotype model, two optional values: 'logistic', 'Gaussian'
#' @param maxiter Numeric, the maximum count for the iterative estimation in when using mixture correlation structure
#' @param tol Numeric, tolerance for the iterative estimation in when using mixture correlation structure
#' @param Hd Numeric, degree of derivative of smooth function which lead to smooth spline of order 2Hd
#' @param VCcorrection Logical variable, indicating whether apply correction for variance components and natural parameters
#' @return This function returns a list object with model parameters and residuals of the NULL model
#' @export


rvmmat_est<-function(y.long, time, y.cov, phe.model = phe.model,maxiter = 50,tol=10^(-6),Hd=2,VCcorrection=FALSE)
{
  
  
  Y<-as.matrix(y.long);
  N<- length(y.long);
  cluster.id<-unique(time[,1]);
  nrow <- m<-length(cluster.id);
  
  time[,2]=time[,2]/max(time[,2])
  
  m=length(unique(time[,1]))
  ncol=ntime=length(unique(time[,2]));  knots=sort(unique(time[,2]))
  Rr=diag(rep(0,ntime))
  for(i in 1:ntime){
    for(jk in i:ntime){
      Rfunc=function(x) (knots[i]-x)^(Hd-1)*(knots[jk]-x)^(Hd-1)/(factorial(Hd-1))^2
      Rr[i,jk]=Rr[jk,i]=stats::integrate(Rfunc,0,knots[i])[[1]]
    }
  }
  Td=matrix(rep(0,Hd*ntime),ncol=Hd)
  for(i in 1:ntime){
    for(jj in 1:Hd){
      Td[i,jj]=(knots[i])^(jj-1)/factorial(jj-1)
    }
  }
  inciN=matrix(rep(0,ntime*N),ncol=ntime)
  for(i in 1:N) inciN[i,which(knots==time[i,2])]=1
  NTd=inciN%*%Td; NRr=inciN%*%Rr; NhalfR=inciN%*%Re(expm::sqrtm(Rr))
  
  
  y.cov=cbind(y.cov,NTd[,-1])
  X <- as.matrix(y.cov);
  X_1 = cbind(1,X)
  
  
  
  fitbi=function(par,vY,mu1){
    solveV=matrix(0,nrow=N,ncol=N)
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    for (i in 1:nrow)
    {
      ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
      mu.i<-mu1[index];
      v1<-matrix(1,length(index), length(index));
      AR.1 <- array(1:ni, dim=c(ni,ni))
      v2 <- 0.7^abs(AR.1-t(AR.1));
      if(ni==1){if(phe.model=="Gaussian"){ v3=par[4];
      } else{  v3<-(mu.i*(1-mu.i))^-1}
        solveV[index,index]<-MASS::ginv(par[1]+par[2]+v3);
      }else{      if(phe.model=="Gaussian"){ v3=diag(par[4], ni);
      } else{  v3<-diag((mu.i*(1-mu.i))^-1)}
        solveV[index,index]<-MASS::ginv(par[1]*v1+par[2]*v2+v3);
      }
      
    }
    for(km in 1:ncol){
      solveVb=solveV%*%NhalfR[,km]
      g=sum(diag(par[3]*t(NhalfR[,km])%*%solveVb))
      solveV=solveV-par[3]*solveVb%*%t(solveVb)/(1+g)
    }
    VX <- solveV%*%X_1
    VY <- solveV%*%vY
    XVX <- crossprod(X_1,VX)
    P=solveV-tcrossprod(VX%*%MASS::ginv(XVX),VX)
    B=MASS::ginv(XVX)%*%crossprod(X_1,VY)
    PY=P%*%vY
    eta=vY-(mu1*(1-mu1))^(-1)*PY
    mu1= inv.logit(eta)
    vY=eta+(mu1*(1-mu1))^(-1)*(y.long-mu1)
    return(list(vY=vY,mu1=mu1,B=B))
  }
  
  
  getD=function(par,vY,mu1,phe.model=phe.model){
    solveV=matrix(0,nrow=N,ncol=N)
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    for (i in 1:nrow)
    {
      ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
      mu.i<-mu1[index];
      v1<-matrix(1,length(index), length(index));
      AR.1 <- array(1:ni, dim=c(ni,ni))
      v2 <- 0.7^abs(AR.1-t(AR.1));
      if(ni==1){if(phe.model=="Gaussian"){ v3=par[4];
      } else{  v3<-(mu.i*(1-mu.i))^-1}
        solveV[index,index]<-MASS::ginv(par[1]+par[2]+v3);
      }else{      if(phe.model=="Gaussian"){ v3=diag(par[4], ni);
      } else{  v3<-diag((mu.i*(1-mu.i))^-1)}
        solveV[index,index]<-MASS::ginv(par[1]*v1+par[2]*v2+v3);
      }
    }
    for(km in 1:ncol){
      solveVb=solveV%*%NhalfR[,km]
      g=sum(diag(par[3]*t(NhalfR[,km])%*%solveVb))
      solveV=solveV-par[3]*solveVb%*%t(solveVb)/(1+g)
    }
    VX <- solveV%*%X_1
    VY <- solveV%*%vY
    XVX <- crossprod(X_1,VX)
    B=MASS::ginv(XVX)%*%crossprod(X_1,VY)
    P=solveV-tcrossprod(VX%*%MASS::ginv(XVX),VX)
    PY=P%*%vY
    if(phe.model=="Gaussian"){
      dql=rep(0,4)
      VPY=matrix(0,ncol=4,nrow=length(PY))
      n.total<-1;n.rep<-as.numeric(table(time[,1]))
      for(i in 1:nrow){
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        AR.1 <- array(1:ni, dim=c(ni,ni))
        R.i<- 0.7^abs(AR.1-t(AR.1))
        dql[1]<-dql[1]+ (sum(PY[index])^2-sum(P[index,index]))/2;
        dql[2]<-dql[2]+(crossprod(PY[index],R.i)%*%PY[index]-sum(diag(P[index,index]%*%R.i)))/2
        VPY[index,1]=rep(sum(PY[index]),ni)
        VPY[index,2]=R.i%*%PY[index]
      }
      dql[3]= (sum(crossprod(PY,NhalfR)^2)-sum(diag(crossprod(NhalfR,P)%*%NhalfR)))/2
      VPY[,3]=NhalfR%*%crossprod(NhalfR,PY)
      dql[4]=(sum(PY^2)-sum(diag(P)))/2
      VPY[,4]=PY
      AI=(crossprod(VPY,P)%*%VPY)/2
      D=MASS::ginv(AI)%*%dql
      
      par1=par+D
      while(any(par1<0)){
        D=D/2
        par1=par+D
        par1[par1<tol&par<tol]=0
      }
      par1[par1<tol]=0
      par=par1
      
    }else{
      dql=rep(0,3)
      VPY=matrix(0,ncol=3,nrow=length(PY))
        n.total<-1;n.rep<-as.numeric(table(time[,1]))
        for(i in 1:nrow){
          ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
          AR.1 <- array(1:ni, dim=c(ni,ni))
          R.i<- 0.7^abs(AR.1-t(AR.1))
          dql[1]<-dql[1]+ (sum(PY[index])^2-sum(P[index,index]))/2;
          dql[2]<-dql[2]+(crossprod(PY[index],R.i)%*%PY[index]-sum(diag(P[index,index]%*%R.i)))/2
          VPY[index,1]=rep(sum(PY[index]),ni)
          VPY[index,2]=R.i%*%PY[index]
        }
      dql[3]= (sum(crossprod(PY,NhalfR)^2)-sum(diag(crossprod(NhalfR,P)%*%NhalfR)))/2
      VPY[,3]=NhalfR%*%crossprod(NhalfR,PY)
      AI=(crossprod(VPY,P)%*%VPY)/2
      D=MASS::ginv(AI)%*%dql
      
      par1=par+D
      while(any(par1<0)){
        D=D/2
        par1=par+D
        par1[par1<tol&par<tol]=0
      }
      par1[par1<tol]=0
      par=par1
      
      fitb=fitbi(par,vY,mu1)
      vY=fitb$vY;mu1=fitb$mu1;B=fitb$B
    }
    return(list(par=par,vY=vY,mu1=mu1,B=B,dql=dql))
  }
  
  
  if(phe.model=="Gaussian"){
    par.init <-  rep(stats::var(y.long, na.rm=T)/4,4);
    par0= par.init
    
    B0= rep(0,dim(X_1)[2]);n=0
    while(n<maxiter){
      dqlAI=getD(par0,y.long,y.long,phe.model)
      parc=dqlAI$par
      B1=dqlAI$B
      #   cat("SIG=",parc, "COV=", B1, "\n");
      
      if(max(abs(c(parc,B1)))>10000) {
        par0=rep(stats::var(y.long, na.rm=T)/4,4)*stats::runif(4);vY0=Y0;mu0=mu
        next
      }
      
      if(max(abs(par0-parc))<=tol*max(abs(par0))&&max(abs(B0-B1))<=tol*max(abs(B0))) break
      n=n+1;B0=B1;par0=parc
    }
    par0=parc;B0=B1;vY0=dqlAI$vY; phi<- par0[4]
  }else{
    fit0 <- stats::glm(Y ~ X, family =  stats::binomial(link = "logit"))
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    Y0 <- eta + fit0$residuals
    mu <- inv.logit(eta)
    par.init <-  rep(stats::var(Y0, na.rm=T)/3,3);
    par0= par.init
    
    B0= rep(0,dim(X_1)[2]);vY0=Y0;mu0=mu;n=0
    while(n<maxiter){
      dqlAI=try(getD(par0,vY0,mu0,phe.model))
      if(class(dqlAI)=="try-error"){
        par0=rep(stats::var(Y0, na.rm=T)/3,3)*stats::runif(3);vY0=Y0;mu0=mu
        next
      }
      parc=dqlAI$par;B1=dqlAI$B;vY0=dqlAI$vY;mu0=dqlAI$mu1;
      
      #  cat("SIG=",parc, "COV=", B1, "\n");
      if(any(is.na(mu0))){
        par0=rep(stats::var(Y0, na.rm=T)/3,3)*stats::runif(3);vY0=Y0;mu0=mu
        next
      }else{ if(min(mu0)==0|max(mu0)==1) {
        par0=rep(stats::var(Y0, na.rm=T)/3,3)*stats::runif(3);vY0=Y0;mu0=mu
        next
      }}
      
      if(max(abs(c(parc,B1)))>10000) {
        par0=rep(stats::var(Y0, na.rm=T)/3,3)*stats::runif(3);vY0=Y0;mu0=mu
        next
      }
      if(max(abs(par0-parc))<=tol*max(abs(par0))&&max(abs(B0-B1))<=tol*max(abs(B0))) break
      n=n+1;B0=B1;par0=parc
    }
    
    par0=dqlAI$par;vY=dqlAI$vY;mu1=dqlAI$mu1; 
    if(max(abs(par0-parc))>tol*max(abs(par0))|max(abs(B0-B1))>tol*max(abs(B0))) warning("Model does not converge!")
    
    if(VCcorrection==TRUE)
    {
      parSS=par0;parSS[1:2]=0;n=0
      while(n<maxiter){
        fitb=fitbi(parSS,vY,mu1)
        B1=fitb$B;vY0=fitb$vY;mu0=fitb$mu1;
        if(max(abs(B0-B1))<tol*max(abs(B0))) break
        n=n+1;B0=B1;vY=vY0;mu1=mu0;
      }
      
      mu.SS=fitb$mu1
      vu = mu.SS*(1-mu.SS);
      dvu = 1 - 2*mu.SS;
      WM0 = vu;
      WM1 = vu*dvu;
      WM2 = -2*vu^2 + (dvu)^2*vu;
      n.total<-1;n.rep<-as.numeric(table(time[,1]))
      ZWZ=Z2WZ2=matrix(0,2,2);XWZ2=matrix(0,ncol(X_1),2)
      for (i in 1:nrow)
      {
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        mu.i<-mu1[index];WM0.i =WM0[index];WM1.i =WM1[index];WM2.i =WM2[index];X.i=X_1[index,]
        v1<-matrix(1,length(index), length(index));
        AR.1 <- array(1:ni, dim=c(ni,ni))
        v2 <- 0.7^abs(AR.1-t(AR.1));
        B=B2=rep(1,ni);
        ihR=expm::sqrtm(v2);ihR2=ihR^2
        ZWZ=ZWZ+matrix(c(sum((t(B)%*%(WM0.i*B))^2),sum((t(B)%*%(WM0.i*ihR))^2),
                         sum((t(ihR)%*%(WM0.i*B))^2),sum((t(ihR)%*%(WM0.i*ihR))^2)),2,2,byrow=T)
        Z2WZ2=Z2WZ2+matrix(c(sum(t(B2)%*%(WM2.i*B2)),sum(t(B2)%*%(WM2.i*ihR2)), 
                             sum(t(ihR2)%*%(WM2.i*B2)),sum(t(ihR2)%*%(WM2.i*ihR2))),2,2,byrow=T)
        if(ni==1){
          XWZ2=XWZ2+cbind(X.i*(WM1.i*B2),X.i*(WM1.i*ihR2[1,1]))
        }else{
          XWZ2=XWZ2+cbind(rowSums(t(X.i)%*%(WM1.i*B2)),rowSums(t(X.i)%*%(WM1.i*ihR2)))
        }
      }
      
      Cp=ZWZ/2
      XWX=t(X_1)%*%(c(WM0)*X_1)
      Cq=Cp+(Z2WZ2-t(XWZ2)%*%MASS::ginv(XWX)%*%XWZ2)/4
      par0[1]=abs(MASS::ginv(Cq[1,1])%*%Cp[1,1]%*%par0[1])
      vY=dqlAI$vY;mu1=dqlAI$mu1;n=0
      while(n<maxiter){
        fitb=fitbi(par0,vY,mu1)
        B1=fitb$B;vY0=fitb$vY;mu0=fitb$mu1;
        if(max(abs(B0-B1))<tol*max(abs(B0))) break
        n=n+1;B0=B1;vY=vY0;mu1=mu0;
      } 
    }
    B0=B1;vY=vY0;mu1=mu0;phi<- 1
  }
  
  
  cat("SIG=",par0, "COV=", B0, "\n");
  tau <- par0
  
  beta<- c(B0)
  mu <- Y.res <- rep(0,length(y.long));
  #V, V.inv, P1, P2
  solveV=matrix(0,nrow=N,ncol=N);V=matrix(0,nrow=N,ncol=N)
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  for (i in 1:nrow){
    ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
    
    v1<-matrix(1,length(index), length(index));
    AR.1 <- array(1:ni, dim=c(ni,ni))
    v2 <- 0.7^abs(AR.1-t(AR.1));
    
    if(ni==1){
      if(phe.model=="Gaussian"){
        v3=tau[4];
      } else{  mu.i<-mu1[index]; v3<-(mu.i*(1-mu.i))^-1
      }
      V[index,index]<- tau[1]+ tau[2]+v3
    }else{
      if(phe.model=="Gaussian"){ v3=diag(tau[4], ni);
      } else{  mu.i<-mu1[index]; v3<-diag((mu.i*(1-mu.i))^-1)
      }
      V[index,index]<- tau[1]*v1+ tau[2]*v2+v3
    }
    solveV[index,index]<-solve(V[index,index]);
  }
  #V=tau[3]* NhalfR%*%t(NhalfR)+V
  for(km in 1:ncol){
    solveVb=solveV%*%NhalfR[,km]
    g=sum(diag(tau[3]*t(NhalfR[,km])%*%solveVb))
    solveV=solveV-tau[3]*solveVb%*%t(solveVb)/(1+g)
  }
  V.inv=solveV
  P1 <- crossprod(X_1,solveV)
  P2<-X_1%*%MASS::ginv(P1%*%X_1)
  if(phe.model=="Gaussian"){
    y.delt <- Y -  X_1 %*% B0
    Y.res=tau[4]*solveV%*%y.delt
    mu <- Y-Y.res
  }else{
    mu=mu1
    Y.res=Y-mu
  }
  
  return(list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,coef = beta, tau = tau,Td=Td,Rr=Rr,inciN=inciN,Hd=Hd,
              est.type = "GLMM", P1=P1, P2=P2,  V.inv = V.inv,phi=phi, Y.res=Y.res, mu = mu,family = stats::binomial()))
}

#'  NULL model estimation for GMMAT and RGMMAT
#'
#' This function estimate the parameters and residuals for the NULL generalized linear mixed effect model
#' @param y.long Long-formatted phenotype vector
#' @param time Time covarites matched with phenotype vector
#' @param y.cov Covariate matrix denoting the covariate variables measured at each time
#' @param phe.model String, the phenotype model, two optional values: 'logistic', 'Gaussian'
#' @param maxiter Numeric, the maximum count for the iterative estimation in when using mixture correlation structure
#' @param tol Numeric, tolerance for the iterative estimation in when using mixture correlation structure
#' @return This function returns a list object with model parameters and residuals of the NULL generalized linear mixed effect model
#' @export

glmm_est<-function(y.long, time, y.cov, phe.model = phe.model,maxiter = 50,tol=10^(-6))
{
  
  
  Y<-as.matrix(y.long);
  N<- length(y.long);
  cluster.id<-unique(time[,1]);
  nrow <- m<-length(cluster.id);
  
  time[,2]=(time[,2]-min(time[,2]))/(max(time[,2])-min(time[,2]))
  
  
  y.cov=cbind(y.cov,time[,2])
  X <- as.matrix(y.cov);
  X_1 = cbind(1,X)
  
  
  
  fitbi=function(par,vY,mu1){
    solveV=list()
    VX=matrix(0,nrow=dim(X_1)[1],ncol=dim(X_1)[2])
    VY=rep(0,length(vY))
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    solveVM=matrix(0,nrow=N,ncol=N)
    for (i in 1:nrow)
    {
      ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
      mu.i<-mu1[index];
      v1<-matrix(1,length(index), length(index));
      AR.1 <- array(1:ni, dim=c(ni,ni))
      v2 <- 0.7^abs(AR.1-t(AR.1));
      if(ni==1){if(phe.model=="Gaussian"){ v3=par[3];
      } else{  v3<-(mu.i*(1-mu.i))^-1}
        solveV[[i]]<-solve(par[1]+par[2]+v3);
        solveVM[index,index]<-solve(par[1]+par[2]+v3);
      }else{      if(phe.model=="Gaussian"){ v3=diag(par[3], ni);
      } else{  v3<-diag((mu.i*(1-mu.i))^-1)}
        solveV[[i]]<-solve(par[1]*v1+par[2]*v2+v3);
        solveVM[index,index]<-solve(par[1]*v1+par[2]*v2+v3);
      }
      VX[index,]=solveV[[i]]%*%X_1[index,]
      VY[index]=solveV[[i]]%*%vY[index]
    }
    
    XVX <- crossprod(X_1,VX)
    P=solveVM-tcrossprod(VX%*%solve(XVX),VX)
    B=solve(XVX)%*%crossprod(X_1,VY)
    PY=P%*%vY
    eta=vY-(mu1*(1-mu1))^(-1)*PY
    mu1= inv.logit(eta)
    vY=eta+(mu1*(1-mu1))^(-1)*(y.long-mu1)
    return(list(vY=vY,mu1=mu1,B=B))
  }
  
  
  getD=function(par,vY,mu1,phe.model=phe.model){
    solveV=list()
    VX=matrix(0,nrow=dim(X_1)[1],ncol=dim(X_1)[2])
    VY=rep(0,length(vY))
    n.total<-1;n.rep<-as.numeric(table(time[,1]))
    solveVM=matrix(0,nrow=N,ncol=N)
    for (i in 1:nrow)
    {
      ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
      mu.i<-mu1[index];
      v1<-matrix(1,length(index), length(index));
      AR.1 <- array(1:ni, dim=c(ni,ni))
      v2 <- 0.7^abs(AR.1-t(AR.1));
      if(ni==1){if(phe.model=="Gaussian"){ v3=par[3];
      } else{  v3<-(mu.i*(1-mu.i))^-1}
        solveV[[i]]<-solve(par[1]+par[2]+v3);
        solveVM[index,index]<-solve(par[1]+par[2]+v3);
      }else{      if(phe.model=="Gaussian"){ v3=diag(par[3], ni);
      } else{  v3<-diag((mu.i*(1-mu.i))^-1)}
        solveV[[i]]<-solve(par[1]*v1+par[2]*v2+v3);
        solveVM[index,index]<-solve(par[1]*v1+par[2]*v2+v3);
      }
      VX[index,]=solveV[[i]]%*%X_1[index,]
      VY[index]=solveV[[i]]%*%vY[index] 
      
    }
    
    XVX <- crossprod(X_1,VX)
    B=solve(XVX)%*%crossprod(X_1,VY)
    P=solveVM-tcrossprod(VX%*%solve(XVX),VX)
    PY=P%*%vY
    if(phe.model=="Gaussian"){ 
      dql=rep(0,3)
      VPY=matrix(0,ncol=3,nrow=length(PY))
      n.total<-1;n.rep<-as.numeric(table(time[,1]))
      for(i in 1:nrow){
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        AR.1 <- array(1:ni, dim=c(ni,ni))
        R.i<- 0.7^abs(AR.1-t(AR.1))
        dql[1]<-dql[1]+ (sum(PY[index])^2-sum(P[index,index]))/2;
        dql[2]<-dql[2]+(crossprod(PY[index],R.i)%*%PY[index]-sum(diag(P[index,index]%*%R.i)))/2
        VPY[index,1]=rep(sum(PY[index]),ni)
        VPY[index,2]=R.i%*%PY[index]
      }
      dql[3]=(sum(PY^2)-sum(diag(P)))/2
      VPY[,3]=PY
      AI=(crossprod(VPY,P)%*%VPY)/2
      D=MASS::ginv(AI)%*%dql
      
      par1=par+D
      while(any(par1<0)){
        D=D/2
        par1=par+D
        par1[par1<tol&par<tol]=0
      }
      par1[par1<tol]=0
      par=par1
      
    }else{
      dql=rep(0,2)
      VPY=matrix(0,ncol=2,nrow=length(PY))
      n.total<-1;n.rep<-as.numeric(table(time[,1]))
      for(i in 1:nrow){
        ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
        AR.1 <- array(1:ni, dim=c(ni,ni))
        R.i<- 0.7^abs(AR.1-t(AR.1))
        dql[1]<-dql[1]+ (sum(PY[index])^2-sum(P[index,index]))/2;
        dql[2]<-dql[2]+(crossprod(PY[index],R.i)%*%PY[index]-sum(diag(P[index,index]%*%R.i)))/2
        VPY[index,1]=rep(sum(PY[index]),ni)
        VPY[index,2]=R.i%*%PY[index]
      }
      AI=(crossprod(VPY,P)%*%VPY)/2
      D=MASS::ginv(AI)%*%dql
      
      par1=par+D
      while(any(par1<0)){
        D=D/2
        par1=par+D
        par1[par1<tol&par<tol]=0
      }
      par1[par1<tol]=0
      par=par1
      
      fitb=fitbi(par,vY,mu1)
      vY=fitb$vY;mu1=fitb$mu1;B=fitb$B
    }
    return(list(par=par,vY=vY,mu1=mu1,B=B,dql=dql))
  }
  
  
  if(phe.model=="Gaussian"){ 
    par.init <-  rep(stats::var(y.long, na.rm=T)/3,3);
    par0= par.init
    
    B0= rep(0,dim(X_1)[2]);n=0
    while(n<maxiter){
      dqlAI=getD(par0,y.long,y.long,phe.model)
      parc=dqlAI$par
      B1=dqlAI$B
      
      if(max(abs(c(parc,B1)))>10000) {
        par0=rep(stats::var(y.long, na.rm=T)/3,3)*stats::runif(3);vY0=Y0;mu0=mu
        next
      }
      
      if(max(abs(par0-parc))<tol*max(abs(par0))&&max(abs(B0-B1))<tol*max(abs(B0))) break
      n=n+1;B0=B1;par0=parc
    } 
    par0=parc;B0=B1;vY0=dqlAI$vY; phi<- par0[3]
  }else{
    fit0 <- stats::glm(Y ~ X, family =  stats::binomial(link = "logit"))
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    Y0 <- eta + fit0$residuals 
    mu <- inv.logit(eta)
    par.init <-  rep(stats::var(Y0, na.rm=T)/2,2);
    par0= par.init
    B0= rep(0,dim(X_1)[2]);vY0=Y0;mu0=mu;n=0
    while(n<maxiter){
      dqlAI=try(getD(par0,vY0,mu0,phe.model))
      if(class(dqlAI)=="try-error"){
        par0=rep(stats::var(Y0, na.rm=T)/2,2)*stats::runif(2);vY0=Y0;mu0=mu
        next
      }
      parc=dqlAI$par;B1=dqlAI$B;vY0=dqlAI$vY;mu0=dqlAI$mu1;
      
      if(any(is.na(mu0))){
        par0=rep(stats::var(Y0, na.rm=T)/2,2)*stats::runif(2);vY0=Y0;mu0=mu
        next
      }else{ if(min(mu0)==0|max(mu0)==1) {
        par0=rep(stats::var(Y0, na.rm=T)/2,2)*stats::runif(2);vY0=Y0;mu0=mu
        next
      }}
      
      if(max(abs(c(parc,B1)))>10000) {
        par0=rep(stats::var(Y0, na.rm=T)/2,2)*stats::runif(2);vY0=Y0;mu0=mu
        next
      }
      if(max(abs(par0-parc))<tol*max(abs(par0))&&max(abs(B0-B1))<tol*max(abs(B0))) break
      n=n+1;B0=B1;par0=parc
    } 
    
    par0=dqlAI$par;vY=dqlAI$vY;mu1=dqlAI$mu1;
    
    B0=B1;vY=vY0;mu1=mu0;
    phi<- 1
  }
  
  cat("SIG=",par0, "COV=", B0, "\n"); 
  tau <- par0
  
  beta<- c(B0)
  mu <- Y.res <- rep(0,length(y.long));
  #V, V.inv, P1, P2
  solveV=matrix(0,nrow=N,ncol=N);V=matrix(0,nrow=N,ncol=N)
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  V=list();V.inv=list()
  solveVM=matrix(0,nrow=N,ncol=N)
  for (i in 1:nrow){
    ni<-n.rep[i];index<-n.total:(n.total+ni-1);n.total<-n.total + ni
    
    v1<-matrix(1,length(index), length(index));
    AR.1 <- array(1:ni, dim=c(ni,ni))
    v2 <- 0.7^abs(AR.1-t(AR.1));
    
    if(ni==1){
      if(phe.model=="Gaussian"){ 
        v3=tau[3];
      } else{  mu.i<-mu1[index]; v3<-(mu.i*(1-mu.i))^-1
      }
      V[[i]]<- tau[1]+ tau[2]+v3
      solveVM[index,index]<-solve(tau[1]+ tau[2]+v3);
    }else{     
      if(phe.model=="Gaussian"){ v3=diag(tau[4], ni);
      } else{  mu.i<-mu1[index]; v3<-diag((mu.i*(1-mu.i))^-1)
      }
      V[[i]]<- tau[1]*v1+ tau[2]*v2+v3
      solveVM[index,index]<-solve(tau[1]*v1+ tau[2]*v2+v3);
    }
    V.inv[[i]]<-solve(V[[i]]);
  }
  
  P1 <- crossprod(X_1,solveVM)
  P2<-X_1%*%MASS::ginv(P1%*%X_1)
  if(phe.model=="Gaussian"){ 
    y.delt <- Y -  X_1 %*% B0
    Y.res=tau[3]*solveVM%*%y.delt
    mu <- Y-Y.res
  }else{
    mu=mu1
    Y.res=Y-mu
  }
  
  return(list(Y=Y,time=time,X=X_1,cluster.id=cluster.id,m=m,coef = beta, tau = tau,est.type = "GLMM", P1=P1, P2=P2, V = V, V.inv = V.inv,phi=phi, Y.res=Y.res, mu = mu,family = binomial()))
}



