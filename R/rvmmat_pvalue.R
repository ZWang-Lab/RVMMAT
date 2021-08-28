#' Calculate prospective and retrospective P-values for RVMMAT test
#'
#' This function tests a SNPs for a given SNP set for a given estimated null model.
#'
#' @param rvmmat.est The output of function "rvmmat_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants.
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#'
#' @return This function returns a dataframe. The row name is the SNP ID, the first column is the prospective score statistics, the second colum is the retrospective score statistics, the third column is the prospective pvalue and the forth column is the restrospective pvalue
#'
#' @export


rvmmat_test <-function(rvmmat.est, G, impute.method='fixed', GRM = NULL)
{
  res<-rvmmat.est$Y.res; phi=rvmmat.est$phi; V.inv<-rvmmat.est$V.inv;X<-rvmmat.est$X;N<-nrow(X)
  m<-rvmmat.est$m;time<-rvmmat.est$time;mu<-rvmmat.est$mu;tau<-rvmmat.est$tau;cluster.id<-rvmmat.est$cluster.id
  snp.names<-colnames(G); family = rvmmat.est$family;Rr=rvmmat.est$Rr;inciN=rvmmat.est$inciN;Td=rvmmat.est$Td;Hd=rvmmat.est$Hd
  ntime=length(unique(time[,2]))
  
  if(!is.null(GRM)){
    dGRM=sqrt(diag(GRM))
    for(i in 1:dim(GRM)[1]) {
      GRM[,i]=GRM[,i]/dGRM[i]
      GRM[i,]=GRM[i,]/dGRM[i]
    }
    
  }
  
  if(is.vector(G))
  {
    G[G==9]<-NA;G<-Impute(G,'fixed')
    center.G<-as.matrix(G-mean(G))
    var_g<-t(center.G)%*%(center.G)/(m-1);
    maf<-min(mean(G)/2, 1-mean(G)/2);
    G<-as.matrix(center.G[match(time[,1],cluster.id),1])
    p=1
  }else
  {
    G[G==9]<-NA
    N_MISS<-sum(is.na(G))
    if(N_MISS>0)
    {
      msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
      warning(msg,call.=F)
      G<-Impute(G,'fixed')
    }
    center.G<-t(t(G)-colMeans(G));
    var_g<-apply(G, 2, stats::var)
    maf<-apply(G, 2, function(x) min(mean(x)/2, 1- mean(x)/2) );
    G<-as.matrix(center.G[match(time[,1],cluster.id),])
    p=dim(G)[2]
  }
  
  P1<-rvmmat.est$P1; P2<-rvmmat.est$P2;
  
  ResN=c(res)*inciN;  ResNTd=ResN%*%Td
  BTResNTd=matrix(0,nrow=m,ncol=Hd)
  
  IDtable=table(time[,1])
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNTd[ss,]= ResNTd[which(time[,1]==ss),]
    }else{ BTResNTd[ss,]= colSums(ResNTd[which(time[,1]==ss),])
    }
  }
  
  halfR=Re(expm::sqrtm(Rr)); ResNR=ResN%*%halfR;
  BTResNR=matrix(rep(0,m*ntime),ncol=ntime)
  
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNR[ss,]= ResNR[which(time[,1]==ss),]
    }else{ BTResNR[ss,]= colSums(ResNR[which(time[,1]==ss),])
    }
  }
  
  BTResNR2colS=colSums(BTResNR^2)
  
  
  if(is.null(GRM)){Var.retro=t(BTResNTd)%*%BTResNTd
  V.Rretro=t(BTResNR)%*%BTResNR
  }else{       Var.retro=t(BTResNTd)%*%GRM%*%BTResNTd
  V.Rretro=t(BTResNR)%*%GRM%*%BTResNR}
  
  type1result=NULL
  
  for( j in 1:p){
    if(is.vector(G)) {
      GN=G*inciN
    }else{
      GN=G[,j]*inciN
    }
    GNTd=GN%*%Td; colnames(GNTd)=1:Hd; GNR=GN%*%Rr
    tranG=GN%*%halfR
    ZGNTd=GNTd- P2%*%(P1%*%GNTd); ZGNR=tranG-P2%*%(P1%*%tranG)
    
    V.pro1<-t(ZGNTd)%*%V.inv%*%ZGNTd;
    
    V.Rr=t(ZGNR)%*%V.inv%*%ZGNR;
    
    tran_res<-res
    
    score2=(t(tran_res)%*%GNTd%*%MASS::ginv(V.pro1)%*%t(GNTd)%*%tran_res)[1,1]/(phi^2)
    
    score2.retro=(t(tran_res)%*%GNTd%*%MASS::ginv(Var.retro*var_g[j])%*%t(GNTd)%*%tran_res)[1,1]
    
    pSmoothsL2=stats::pchisq(score2,df=Hd,lower.tail=F)
    pSmoothsL.retro2=stats::pchisq(score2.retro,df=Hd,lower.tail=F)
    
    HalfTSKAT=t(tranG)%*%tran_res
    TSKAT=sum(HalfTSKAT^2)
    lambdaS=eigen(V.Rr,symmetric = TRUE,only.values = TRUE)$values
    
    if(sum(lambdaS)==0) { pSmoothsSKAT2=1
    }else{ pSmoothsSKAT2=generalchisq(lambdaS,TSKAT/(phi^2))}
    
    lambdaS.retro=eigen(V.Rretro*var_g[j],symmetric = TRUE,only.values = TRUE)$values
    if(sum( lambdaS.retro)==0) {pSmoothsSKAT.retro2=1
    }else{
      pSmoothsSKAT.retro2=generalchisq(lambdaS.retro,TSKAT)}
    
    Tp=(tan((0.5-pSmoothsL2)*pi)+tan((0.5-pSmoothsSKAT2)*pi))/2
    pSmoothsACAT2=0.5-atan(Tp)/pi
    Tr=(tan((0.5-pSmoothsL.retro2)*pi)+tan((0.5-pSmoothsSKAT.retro2)*pi))/2
    pSmoothsACAT.retro2=0.5-atan(Tr)/pi
    
    result=cbind(pSmoothsL2, pSmoothsACAT2, pSmoothsL.retro2,  pSmoothsACAT.retro2)
    type1result=rbind(type1result,result)
  }
  return(type1result)
  
}





Impute<-function(Z, impute.method){
  if(is.vector(Z)){
    if(impute.method =="random"){
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-stats::rbinom(length(IDX),2,maf1)
      }
    } else if(impute.method =="fixed"){
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-2 * maf1
      }
    } else if(impute.method =="bestguess") {
      IDX<-which(is.na(Z))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX])/2
        Z[IDX]<-round(2 * maf1)
      }
    } else {
      stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
    }
  } else{
    p<-dim(Z)[2]
    if(impute.method =="random"){
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-stats::rbinom(length(IDX),2,maf1)
        }
      }
    } else if(impute.method =="fixed"){
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-2 * maf1
        }
      }
    } else if(impute.method =="bestguess") {
      for(i in 1:p){
        IDX<-which(is.na(Z[,i]))
        if(length(IDX) > 0){
          maf1<-mean(Z[-IDX,i])/2
          Z[IDX,i]<-round(2 * maf1)
        }
      }
    } else {
      stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
    }
  }
  return(as.matrix(Z))
}




generalchisq=function(lambda,Q){
  muq=sum(lambda)
  sigmaq=sqrt(2*sum(lambda^2))
  s1=sum(lambda^3)/sum(lambda^2)^1.5
  s2=sum(lambda^4)/sum(lambda^2)^2
  if((s1^2)>s2){
    a=1/(s1-sqrt(s1^2-s2))
    delta=s1*a^3-a^2
    l=a^2-2*delta
  }else{
    delta=0
    l=sum(lambda^2)^3/sum(lambda^3)^2
  }
  stats::pchisq((sum(Q)-muq)/sigmaq*sqrt(2*(l+2*delta))+l+delta,df=l, ncp=delta,lower.tail=FALSE)}

