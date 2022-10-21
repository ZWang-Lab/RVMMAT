library(data.table)

type1c=function(vect,x){
ln=length(x)
result=rep(0,ln)
 for(i in 1:ln){
  result[i]=sum(vect<x[i])
  }
return(result)
}

powerc=function(vect,x){
sum(vect<x)
}



file.path="RVMMAT_simu\\"


#type1 logistic
logtype1r=fread(paste0(file.path,"logistic\\type1\\random\\type1result.txt"),header=TRUE)
logtype1b=fread(paste0(file.path,"logistic\\type1\\baseline\\type1result.txt"),header=TRUE)
logtype1s=fread(paste0(file.path,"logistic\\type1\\sum\\type1result.txt"),header=TRUE)
#liability
liatype1r=fread(paste0(file.path,"liability\\type1\\random\\type1result.txt"),header=TRUE)
liatype1b=fread(paste0(file.path,"liability\\type1\\baseline\\type1result.txt"),header=TRUE)
liatype1s=fread(paste0(file.path,"liability\\type1\\sum\\type1result.txt"),header=TRUE)

dim(logtype1r)
dim(logtype1b)
dim(logtype1s)
dim(liatype1r)
dim(liatype1b)
dim(liatype1s)

plotn=dim(logtype1r)[1]
alpha=c(0.01,0.001,0.0001)

type1=(cbind(rep(alpha*plotn,2),c(apply(logtype1r[,],2,type1c,x=alpha)),c(apply(logtype1b[,],2,type1c,x=alpha)),c(apply(logtype1s[,],2,type1c,x=alpha)),
c(apply(liatype1r[,],2,type1c,x=alpha)),c(apply(liatype1b[,],2,type1c,x=alpha)),c(apply(liatype1s[,],2,type1c,x=alpha)))/plotn)

type1=type1[c(4:6,1:3),]
colnames(type1)=c("Nominal Level",rep(c("random","baseline","sum"),2))
methods=c("RVMMAT","VMMAT")
rownames(type1)=rep(methods,each=3)

type1

c(0.01-1.96*sqrt(0.01*0.99/plotn),0.01+1.96*sqrt(0.01*0.99/plotn))
c(0.001-1.96*sqrt(0.001*0.999/plotn),0.001+1.96*sqrt(0.001*0.999/plotn))
c(0.0001-1.96*sqrt(0.0001*0.9999/plotn),0.0001+1.96*sqrt(0.0001*0.9999/plotn))


#power logistic (9-24)
logpowerb=read.table(paste0(file.path,"logistic\\power\\baseline\\type1result.txt"),header=TRUE)
logpowerb2=read.table(paste0(file.path,"logistic\\power\\baseline\\type1result2.txt"),header=TRUE)
logpowerb3=read.table(paste0(file.path,"logistic\\power\\baseline\\type1result3.txt"),header=TRUE)
logpowerb4=read.table(paste0(file.path,"logistic\\power\\baseline\\type1result4.txt"),header=TRUE)

logpowerr=read.table(paste0(file.path,"logistic\\power\\random\\type1result.txt"),header=TRUE)
logpowerr2=read.table(paste0(file.path,"logistic\\power\\random\\type1result2.txt"),header=TRUE)
logpowerr3=read.table(paste0(file.path,"logistic\\power\\random\\type1result3.txt"),header=TRUE)
logpowerr4=read.table(paste0(file.path,"logistic\\power\\random\\type1result4.txt"),header=TRUE)

logpowers=read.table(paste0(file.path,"logistic\\power\\sum\\type1result.txt"),header=TRUE)
logpowers2=read.table(paste0(file.path,"logistic\\power\\sum\\type1result2.txt"),header=TRUE)
logpowers3=read.table(paste0(file.path,"logistic\\power\\sum\\type1result3.txt"),header=TRUE)
logpowers4=read.table(paste0(file.path,"logistic\\power\\sum\\type1result4.txt"),header=TRUE)

#power liability
liapowerb=read.table(paste0(file.path,"liability\\power\\baseline\\type1result.txt"),header=TRUE)
liapowerb2=read.table(paste0(file.path,"liability\\power\\baseline\\type1result2.txt"),header=TRUE)
liapowerb3=read.table(paste0(file.path,"liability\\power\\baseline\\type1result3.txt"),header=TRUE)
liapowerb4=read.table(paste0(file.path,"liability\\power\\baseline\\type1result4.txt"),header=TRUE)

liapowerr=read.table(paste0(file.path,"liability\\power\\random\\type1result.txt"),header=TRUE)
liapowerr2=read.table(paste0(file.path,"liability\\power\\random\\type1result2.txt"),header=TRUE)
liapowerr3=read.table(paste0(file.path,"liability\\power\\random\\type1result3.txt"),header=TRUE)
liapowerr4=read.table(paste0(file.path,"liability\\power\\random\\type1result4.txt"),header=TRUE)

liapowers=read.table(paste0(file.path,"liability\\power\\sum\\type1result.txt"),header=TRUE)
liapowers2=read.table(paste0(file.path,"liability\\power\\sum\\type1result2.txt"),header=TRUE)
liapowers3=read.table(paste0(file.path,"liability\\power\\sum\\type1result3.txt"),header=TRUE)
liapowers4=read.table(paste0(file.path,"liability\\power\\sum\\type1result4.txt"),header=TRUE)


dim(liapowerb);dim(liapowerb2);dim(liapowerb3);dim(liapowerb4)
dim(liapowerr);dim(liapowerr2);dim(liapowerr3);dim(liapowerr4)
dim(liapowers);dim(liapowers2);dim(liapowers3);dim(liapowers4)


dim(logpowerb);dim(logpowerb2);dim(logpowerb3);dim(logpowerb4)
dim(logpowerr);dim(logpowerr2);dim(logpowerr3);dim(logpowerr4)
dim(logpowers);dim(logpowers2);dim(logpowers3);dim(logpowers4)

plotn=dim(liapowerb2)[1]
alpha=10^-3

#
logisb=(rbind(apply(logpowerb,2,powerc,x=alpha),apply(logpowerb2,2,powerc,x=alpha),
apply(logpowerb3,2,powerc,x=alpha),apply(logpowerb4,2,powerc,x=alpha))/plotn)

logisr=(rbind(apply(logpowerr,2,powerc,x=alpha),apply(logpowerr2,2,powerc,x=alpha),
apply(logpowerr3,2,powerc,x=alpha),apply(logpowerr4,2,powerc,x=alpha))/plotn)

logiss=(rbind(apply(logpowers,2,powerc,x=alpha),apply(logpowers2,2,powerc,x=alpha),
apply(logpowers3,2,powerc,x=alpha),apply(logpowers4,2,powerc,x=alpha))/plotn)


liab=(rbind(apply(liapowerb,2,powerc,x=alpha),apply(liapowerb2,2,powerc,x=alpha),
apply(liapowerb3,2,powerc,x=alpha),apply(liapowerb4,2,powerc,x=alpha))/plotn)


liar=(rbind(apply(liapowerr,2,powerc,x=alpha),apply(liapowerr2,2,powerc,x=alpha),
apply(liapowerr3,2,powerc,x=alpha),apply(liapowerr4,2,powerc,x=alpha))/plotn)

lias=(rbind(apply(liapowers,2,powerc,x=alpha),apply(liapowers2,2,powerc,x=alpha),
apply(liapowers3,2,powerc,x=alpha),apply(liapowers4,2,powerc,x=alpha))/plotn)
#rownames(logis)=c("random","baseline","sum")


xname5=c(0.6,0.63,0.66,0.69)
colorm=c("blue","blue","red","red","green")
pchm=c(2,1,2,1,2)


tiff(paste0(file.path,"PowerVCC.tiff"), width = 12, height = 8,compression = "lzw", units = 'in', res = 300)

par(mfrow=c(2,3),mar=c(5,5,4,1),oma=c(4,0,3,0),cex.axis=1.9,cex.lab=1.9,cex.main=2)
 matplot(xname5,logisr, main = "Random",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2)
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)
title(ylab=expression("Power"), line=3, cex.lab=1.9)

 matplot(xname5,logisb, main = "Baseline",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2)
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)

mtext(expression(bold("Liability Threshold Model")),side = 1,line=6.5,outer=FALSE,cex=2)
mtext(expression(bold("Logistic Mixed Model")),outer=TRUE,line=0,cex=2)

 matplot(xname5,logiss, main = "Sum",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2)
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)

 matplot(xname5,liar, main = "",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2) 
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)
title(ylab=expression("Power"), line=3, cex.lab=1.9)
title(xlab=expression(gamma), line=4.1, cex.lab=2.5)

 matplot(xname5,liab, main = "",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2) 
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)
title(xlab=expression(gamma), line=4.1, cex.lab=2.5)

 matplot(xname5,lias, main = "",xaxt="n",ylim=c(0,0.9),
ylab="",xlab="",type = c("l"),col = colorm,lty=pchm,bty="l",lwd=2)
axis(side=1,at=c(0.5,0.7),lwd=2)
axis(side=1,at=xname5)
axis(2,at=c(-0.1,1),lwd=2)
title(xlab=expression(gamma), line=4.1, cex.lab=2.5)

methods=c("RVMMAT","VMMAT","Copula","RGMMAT","GMMAT")
colorm2=c("red","red","Green","blue","blue")
pchm2=c(1,2,2,1,2)

legend(x=0.37,y=-0.3,xpd=NA,x.intersp=1,text.width=0.04, legend=methods,cex=1.9,pt.cex=2,
ncol=5,col=colorm2,bty="n",lty=pchm2,lwd=2)

dev.off()


