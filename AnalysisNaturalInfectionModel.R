###################################
### Analysis Natural Infection
### By JM Prada

## Load data and set working directory
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
rm(list=ls())
load("MainResults.RData")

## Sample from posteriors
set.seed(1)
ndraws <- 1000
IDs <- sample.int(20000,ndraws)

PrevalenceFarm <- array(NA,dim=c(nFarm,nParasite,ndraws))
FKF <- array(NA,dim=c(nFarm,nParasite,nPool,nSamplesPool,ndraws))
Sed <- array(NA,dim=c(nFarm,nParasite,nPool,nSamplesPool,ndraws))

for (f in 1:nFarm){
  for (s in 1:nParasite){
    PrevalenceFarm[f,s,] <- c(MainModel$mcmc[[1]][,paste("P[",f,",",s,"]",sep = "")],
                             MainModel$mcmc[[2]][,paste("P[",f,",",s,"]",sep = "")])[IDs]
    for (p in 1:nPool){
      for(i in 1:nSamplesPool){
        FKF[f,s,p,i,] <- c(MainModel$mcmc[[1]][,paste("FKF[",f,",",s,",",p,",",i,"]",sep = "")],
                          MainModel$mcmc[[2]][,paste("FKF[",f,",",s,",",p,",",i,"]",sep = "")])[IDs]
        Sed[f,s,p,i,] <- c(MainModel$mcmc[[1]][,paste("Sed[",f,",",s,",",p,",",i,"]",sep = "")],
                           MainModel$mcmc[[2]][,paste("Sed[",f,",",s,",",p,",",i,"]",sep = "")])[IDs]
      }
    }
  }
}

## Simulated prevalence
## Fas = Fasciola, Cal = Calicophoron
PrevFas <- apply(PrevalenceFarm[,1,],1,quantile,c(.025,.5,.975))
PrevCal <- apply(PrevalenceFarm[,2,],1,quantile,c(.025,.5,.975))

## Prevalence by Mini-FLOTAC (data)
PrevFasFT <-apply(Flotac[,1,,],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)
PrevCalFT <-apply(Flotac[,2,,],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)

## Prevalence by FKF (estimated by model)
PrevFasFF <-sapply(1:ndraws,function(x) {apply(FKF[,1,,,x],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)})
PrevFasSd <-sapply(1:ndraws,function(x) {apply(Sed[,1,,,x],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)})

## Prevalence by Sed (estimated by model)
PrevCalFF <-sapply(1:ndraws,function(x) {apply(FKF[,2,,,x],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)})
PrevCalSd <-sapply(1:ndraws,function(x) {apply(Sed[,2,,,x],1,function(x)length(which(x!=0)))/(nPool*nSamplesPool)})

## Calculate quantiles for the model estimates (FKF and Sed)
PrevFasFFQuant <- apply(PrevFasFF,1,quantile,c(.025,.5,.975))
PrevFasSdQuant <- apply(PrevFasSd,1,quantile,c(.025,.5,.975))

PrevCalFFQuant <- apply(PrevCalFF,1,quantile,c(.025,.5,.975))
PrevCalSdQuant <- apply(PrevCalSd,1,quantile,c(.025,.5,.975))

spt = .2 # separation between points

pdf("Figure 4 PrevalenceEstimation.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(1,nFarm*nParasite),ylim=c(0,1), axes = F,xlab = "", ylab = "")

arrows(x0=(1:20)+(spt/2+spt),y0=c(PrevFas[1,],PrevCal[1,]),
       y1=c(PrevFas[3,],PrevCal[3,]),angle = 90, code = 3, length = .055)
points((1:20)+(spt/2+spt),c(PrevFas[2,],PrevCal[2,]), pch=19)
points((1:20)-(spt/2+spt),c(PrevFasFT,PrevCalFT), pch=17, col="blueviolet")
arrows(x0=(1:20)-(spt/2),y0=c(PrevFasFFQuant[1,],PrevCalFFQuant[1,]), col = "darkred",
       y1=c(PrevFasFFQuant[3,],PrevCalFFQuant[3,]),angle = 90, code = 3, length = .055)
points((1:20)-(spt/2),c(PrevFasFFQuant[2,],PrevCalFFQuant[2,]), pch=19, col = "darkred")
arrows(x0=(1:20)+(spt/2),y0=c(PrevFasSdQuant[1,],PrevCalSdQuant[1,]), col = "darkgreen",
       y1=c(PrevFasSdQuant[3,],PrevCalSdQuant[3,]),angle = 90, code = 3, length = .055)
points((1:20)+(spt/2),c(PrevFasSdQuant[2,],PrevCalSdQuant[2,]), pch=19,col = "darkgreen")

abline(v=10.5,lty="dashed")

axis(1, at=1:20, labels=rep(1:10,2))
axis(2, at=seq(0,1,by=.2), labels = seq(0,1,by=.2)*100)

mtext("Prevalence (%)",side=2,cex=1,line=1.2)
mtext("Farm number",side=1,cex=1,line=1.2)

text(5.25,1.05,labels = expression(italic("F. Hepatica")), xpd=NA)
text(15.75,1.05,labels = expression(italic("C. daubneyi")), xpd=NA)

legend("topleft",c("Mini-FLOTAC","Flukefinder (sim)","Sedimentation (sim)","Model Estimate"),
       col=c("blueviolet","darkred","darkgreen","black"),
       bty='n',cex=0.75,lty=c(NA,1,1,1),pch = c(17,19,19,19))

dev.off()



#### Plot prevalence distributions for all farms ####
do.plot <- function(f,s){
  Prev <- PrevalenceFarm[f,s,]
  plot(NA,xlim = c(0,1),ylim = c(0,1),axes = F, xlab = "", ylab = "",
       main = paste("Farm ",f,", ", Species[s], sep=""))
  lines(density(Prev)$x,density(Prev)$y/max(density(Prev)$y))
  
  axis(1, at = seq(0,1,by=0.2), labels = seq(0,1,by=0.2)*100)
  axis(2)
  
  mtext("Scaled Density",side=2,cex=1,line=1.2)
  mtext("Prevalence (%)",side=1,cex=1,line=1.2)
}

##
pdf("Prevalence F. hepatica.pdf",width=7*3,height=4.64*4)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(4,3))

for (f in 1:10){
  do.plot(f,1)
}

dev.off()

##
pdf("Prevalence C. daubneyi.pdf",width=7*3,height=4.64*4)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfrow=c(4,3))

for (f in 1:10){
  do.plot(f,2)
}

dev.off()