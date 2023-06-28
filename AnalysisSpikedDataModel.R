##################################
### Analysis Spiked Data Model
### By JM Prada

## Load data and set working directory
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
rm(list=ls())

load("SpikedResults.RData")

## Generate sequence for x axis
nEggs <- seq(1,100,by=.5)*10 #EPG * 10 to make 10 grams

## Sample from posteriors and simulate counts
nRepeats = 1000
set.seed(101)
IDs <- sample.int(20000,nRepeats)
rt_Ft <- c(NoiseSpiked$mcmc[[1]][,"rt_Ft"],NoiseSpiked$mcmc[[2]][,"rt_Ft"])[IDs]
rt_FF <- c(NoiseSpiked$mcmc[[1]][,"rt_FF"],NoiseSpiked$mcmc[[2]][,"rt_FF"])[IDs]
rt_Sd <- c(NoiseSpiked$mcmc[[1]][,"rt_Sd"],NoiseSpiked$mcmc[[2]][,"rt_Sd"])[IDs]

FtCounts <- sapply(1:nRepeats, function(x){rnbinom(length(nEggs), size=rt_Ft[x], mu=nEggs/50)})*5
FFCounts <- sapply(1:nRepeats, function(x){rnbinom(length(nEggs), size=rt_FF[x], mu=nEggs/5)})/2
SdCounts <- sapply(1:nRepeats, function(x){rnbinom(length(nEggs), size=rt_Sd[x], mu=nEggs)})/10

## Calculate quantiles for counts
FtQuant <- apply(FtCounts,1,quantile,c(.025,.5,.975))
FFQuant <- apply(FFCounts,1,quantile,c(.025,.5,.975))
SdQuant <- apply(SdCounts,1,quantile,c(.025,.5,.975))

library(scales)

#### Using quantiles plot Real vs Simulated EPGs ####
pdf("Figure3-SimulatedvsRealEPGs.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(0,max(nEggs/10)),ylim=c(0,max(nEggs/10)*5), axes = F,xlab = "", ylab = "")

polygon(c(nEggs/10,rev(nEggs/10)),c(SdQuant[1,],rev(SdQuant[3,])), col=alpha("darkgreen",.3), border=F)
polygon(c(nEggs/10,rev(nEggs/10)),c(FFQuant[1,],rev(FFQuant[3,])), col=alpha("darkred",.3), border=F)
polygon(c(nEggs/10,rev(nEggs/10)),c(FtQuant[1,],rev(FtQuant[3,])), col=alpha("blueviolet",.3), border=F)

lines(nEggs/10,SdQuant[2,], col="darkgreen",lwd=4)
lines(nEggs/10,FFQuant[2,], col="darkred",lwd=4)
lines(nEggs/10,FtQuant[2,], col="blueviolet",lwd=4)

abline(a=0,b=1,lty="dashed")

axis(1)
axis(2)

mtext("Estimated infection intensity (EPG)",side=2,cex=1,line=1.2)
mtext("True infection intensity (EPG)",side=1,cex=1,line=1.2)

legend("topleft",c("Mini-FLOTAC","Flukefinder","Sedimentation"),
       col=c("blueviolet","darkred","darkgreen"),
       bty='n',cex=0.75,lty=1)

dev.off()

## Calculate quantiles for accuracy (i.e. absolute difference)
FtQuantAc <- apply((abs(FtCounts-nEggs/10)),1,quantile,c(.025,.5,.975))
FFQuantAc <- apply((abs(FFCounts-nEggs/10)),1,quantile,c(.025,.5,.975))
SdQuantAc <- apply((abs(SdCounts-nEggs/10)),1,quantile,c(.025,.5,.975))

#### Plotting as accuracy in the EPG ####
pdf("AccuracyEPGs.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(0,max(nEggs/10)),ylim=c(0,max(nEggs/10)*5), axes = F,xlab = "", ylab = "")

polygon(c(nEggs/10,rev(nEggs/10)),c(SdQuantAc[1,],rev(SdQuantAc[3,])), col=alpha("darkgreen",.3), border=F)
polygon(c(nEggs/10,rev(nEggs/10)),c(FFQuantAc[1,],rev(FFQuantAc[3,])), col=alpha("darkred",.3), border=F)
polygon(c(nEggs/10,rev(nEggs/10)),c(FtQuantAc[1,],rev(FtQuantAc[3,])), col=alpha("blueviolet",.3), border=F)

lines(nEggs/10,SdQuantAc[2,], col="darkgreen",lwd=4)
lines(nEggs/10,FFQuantAc[2,], col="darkred",lwd=4)
lines(nEggs/10,FtQuantAc[2,], col="blueviolet",lwd=4)

abline(h=0,lty="dashed")

axis(1)
axis(2)

mtext("Error in EPG",side=2,cex=1,line=1.2)
mtext("True EPG",side=1,cex=1,line=1.2)

legend("topleft",c("Mini-FLOTAC","Flukefinder","Sedimentation"),
       col=c("blueviolet","darkred","darkgreen"),
       bty='n',cex=0.75,lty=1)

dev.off()


## Calculate expected zero counts
FtPropNeg <- sapply(1:nrow(FtCounts),function(x){length(which(FtCounts[x,]==0))})/nRepeats
FFPropNeg <- sapply(1:nrow(FFCounts),function(x){length(which(FFCounts[x,]==0))})/nRepeats
SdPropNeg <- sapply(1:nrow(SdCounts),function(x){length(which(SdCounts[x,]==0))})/nRepeats

#### Check values real ones ####
#dt <- read.csv("SpikedInfection.csv")
#Dosage <- unique(dt$Dose)
#TFtPropNeg <- sapply(Dosage,function(x){length(which(dt$Mini.FLOTAC...EPG[which(dt$Dose==x)]==0))})/12
#TFFPropNeg <- sapply(Dosage,function(x){length(which(dt$Flukefinder..EPG[which(dt$Dose==x)]==0))})/12
#TSdPropNeg <- sapply(Dosage,function(x){length(which(dt$Sedimentation...EPG[which(dt$Dose==x)]==0))})/12
#####

#### Plot expected zero (negative) counts ####
pdf("ExpectedNegativeCounts.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(0,max(nEggs/10)),ylim=c(0,.3), axes = F,xlab = "", ylab = "")

lines(nEggs/10,FtPropNeg, col="blueviolet",lwd=2)
lines(nEggs/10,SdPropNeg, col="darkgreen",lwd=2)
lines(nEggs/10,FFPropNeg, col="darkred",lwd=2)

axis(1)
axis(2, at=seq(0,.3,by=.05), labels = seq(0,.3,by=.05)*100)

mtext("Expected Negative proportion (%)",side=2,cex=1,line=1.2)
mtext("True EPG",side=1,cex=1,line=1.2)

legend("topright",c("Mini-FLOTAC","Flukefinder","Sedimentation"),
       col=c("blueviolet","darkred","darkgreen"),
       bty='n',cex=0.75,lty=1)

dev.off()

#### Plot Sensitivity ####
pdf("Figure 2 Sensitivity.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(0,max(nEggs/10)),ylim=c(0,1), axes = F,xlab = "", ylab = "")

lines(nEggs/10,1-FtPropNeg, col="blueviolet",lwd=2)
lines(nEggs/10,1-SdPropNeg, col="darkgreen",lwd=2)
lines(nEggs/10,1-FFPropNeg, col="darkred",lwd=2)

#points(Dosage,1-TFtPropNeg, pch=19, col="blueviolet")
#points(Dosage,1-TSdPropNeg, pch=19, col="darkgreen")
#points(Dosage,1-TFFPropNeg, pch=19, col="darkred")

axis(1)
axis(2, at=seq(0,1,by=.2), labels = seq(0,1,by=.2)*100)

mtext("Sensitivity (%)",side=2,cex=1,line=1.2)
mtext("Infection intensity (EPG)",side=1,cex=1,line=1.2)

legend("bottomright",c("Mini-FLOTAC","Flukefinder","Sedimentation"),
       col=c("blueviolet","darkred","darkgreen"),
       bty='n',cex=0.75,lty=1)

dev.off()

######


#### Plot Points ####
pdf("SimulatedvsREalEPGsPoints.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim=c(0,1000),ylim=c(0,5000), axes = F,xlab = "", ylab = "")

points(rep(nEggs/10,nRepeats),SdCounts, col=alpha("darkgreen",.03),pch=19)
points(rep(nEggs/10,nRepeats),FFCounts, col=alpha("darkred",.03),pch=19)
points(rep(nEggs/10,nRepeats),FtCounts, col=alpha("grey",.04),pch=19)

abline(a=0,b=1,lty="dashed")

axis(1)
axis(2)

mtext("Estimated EPG",side=2,cex=1,line=1.2)
mtext("True EPG",side=1,cex=1,line=1.2)

legend("topleft",c("Mini-Flotac","Fluke-Finder","Sedimentation"),
       col=c("grey","darkred","darkgreen"),
       bty='n',cex=0.75,pch=19)

dev.off()
