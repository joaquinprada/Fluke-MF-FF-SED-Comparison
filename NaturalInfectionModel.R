####################################
### Run model for Natural Infection
### By JM Prada

## Load data and set working directory
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
rm(list=ls())
load("SpikedResults.RData")

dtIndi <- read.csv("FieldData.csv")
dtPool <- read.csv("PooledSamples.csv")

## Load posteriors from spiked data
library(fitdistrplus)

## Prior for rt_Ft
fit <- fitdist(c(NoiseSpiked$mcmc[[1]][,"rt_Ft"],NoiseSpiked$mcmc[[2]][,"rt_Ft"]),
               distr = "gamma", method = "mle")
rt_FTparams <- fit$estimate

## Prior for rt_FF
fit <- fitdist(c(NoiseSpiked$mcmc[[1]][,"rt_FF"],NoiseSpiked$mcmc[[2]][,"rt_FF"]),
               distr = "gamma", method = "mle")
rt_FFparams <- fit$estimate

## Prior for rt_Sd
fit <- fitdist(c(NoiseSpiked$mcmc[[1]][,"rt_Sd"],NoiseSpiked$mcmc[[2]][,"rt_Sd"]),
               distr = "gamma", method = "mle")
rt_Sdparams <- fit$estimate


## Data preparation ##
dtIndi$Parasite[which(dtIndi$Parasite=="F. hepatica ")] <- "F. hepatica"
dtIndi$Parasite[which(dtIndi$Parasite=="C. daubneyi ")] <- "C. daubneyi"
Farm <- unique(dtIndi$Farm)
Species <- as.factor(unique(dtIndi$Parasite))
Pool <- unique(dtIndi$Pool)

nFarm <- length(Farm)
nParasite <- length(Species)
nPool <- length(Pool)
nSamplesPool <- nrow(dtIndi)/(nFarm*nParasite*length(Pool))

Flotac <- array(NA,dim=c(nFarm,nParasite,nPool,nSamplesPool))
FlotacPool <- array(NA,dim=c(nFarm,nParasite,nPool))
FKFPool <- array(NA,dim=c(nFarm,nParasite,nPool))
SedPool <- array(NA,dim=c(nFarm,nParasite,nPool))

for( f in Farm){
  for (s in as.numeric(Species)){
    for (p in Pool){
      Flotac[f,s,p,] <- dtIndi$Flotac[which(dtIndi$Farm==f &
                                              dtIndi$Parasite==Species[s] &
                                              dtIndi$Pool==p)]/5 #0.2g sample
      FlotacPool[f,s,p] <- dtPool$Flotac[which(dtPool$Farm==f &
                                                 dtPool$Parasite==Species[s] &
                                                 dtPool$Pool==p)]/5
      FKFPool[f,s,p] <- dtPool$FlukeFinder[which(dtPool$Farm==f &
                                                   dtPool$Parasite==Species[s] &
                                                   dtPool$Pool==p)]*2
      SedPool[f,s,p] <- dtPool$Sedimentation[which(dtPool$Farm==f &
                                                     dtPool$Parasite==Species[s] &
                                                     dtPool$Pool==p)]*10
    }
  }
}

## Provide a distribution for Fasciola as a prior, as it is mostly never there
## Scale to 10g (as that is the biggest sample)
## testing based on data fit <- fitdist(c(Flotac[,,1,])*50, distr = "gamma", method = "mme")
## Mean EPG taken from Rinaldi et al. Geospatial Health 9(2), 2015, pp. 309-317
## mean=27.8 25-50-75th quantiles (10.0 - 21.0 - 48.0)
refQ <- c(10,21,48)*10 #scale EPGs to 10grams
fitGamma <- function(x) sum(abs(refQ-qgamma(c(.25,.5,.75),shape=x[1],rate=x[2]))^2)
fit <- optim(c(.5,.1),fitGamma)
PriorFasciola <- fit$par

## Set seed ##
.RNG.seed <- function(chain)
  return( switch(chain, "1"= 1, "2"= 2) )
.RNG.name <- function(chain)
  return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )

## Set initial status to 1
status <- array(1,dim=c(nFarm,nParasite,nPool,nSamplesPool))

## Define model
m <- "model {
  for (f in 1:nFarm){
    for (s in 1:nParasite){
      
      P[f,s] ~ dbeta(1,1)
      
      for (p in 1:nPool){
        for(i in 1:nSamplesPool){
          status[f,s,p,i] ~ dbern(P[f,s])
        
          nEggs[f,s,p,i] ~ dgamma(sh[s],rt[s])
          
          Flotac[f,s,p,i] ~ dnegbin( (rt_Ft/(nEggs[f,s,p,i]*status[f,s,p,i]/50+rt_Ft)), rt_Ft)
          FKF[f,s,p,i] ~ dnegbin( (rt_FF/(nEggs[f,s,p,i]*status[f,s,p,i]/5+rt_FF)), rt_FF)
          Sed[f,s,p,i] ~ dnegbin( (rt_Sd/(nEggs[f,s,p,i]*status[f,s,p,i]+rt_Sd)), rt_Sd)
        }
        
        nEggsPool[f,s,p] <- sum(nEggs[f,s,p,]*status[f,s,p,])/nSamplesPool
      
        ## Flotac
        FlotacPool[f,s,p] ~ dnegbin( (rt_Ft/(nEggsPool[f,s,p]/50+rt_Ft)), rt_Ft)
    
        ## Flukefinder
        FKFPool[f,s,p] ~ dnegbin( (rt_FF/(nEggsPool[f,s,p]/5+rt_FF)), rt_FF)
    
        ## Sedimentation
        SedPool[f,s,p] ~ dnegbin( (rt_Sd/(nEggsPool[f,s,p]+rt_Sd)), rt_Sd)
      }
    }
  }
  
  ## Priors at parasite level
  sh[1] <- PriorFasciola[1]
  rt[1] <- PriorFasciola[2]
  sh[2] ~ dgamma(0.001,0.001)
  mn[2] ~ dgamma(0.001,0.001)
  rt[2] <- sh[2]/mn[2]

  rt_Ft ~ dgamma(rt_FTparams[1],rt_FTparams[2])
  rt_FF ~ dgamma(rt_FFparams[1],rt_FFparams[2])
  rt_Sd ~ dgamma(rt_Sdparams[1],rt_Sdparams[2])
  
  #inits# .RNG.seed, .RNG.name, status
  #data# nFarm, nParasite, nPool, nSamplesPool, Flotac, FlotacPool, FKFPool, SedPool, rt_FTparams, rt_FFparams, rt_Sdparams, PriorFasciola
  #monitor# rt_Ft, rt_FF, rt_Sd, P, sh, rt, FKF, Sed
}"


MainModel <- run.jags(m, burnin=1000, sample=10000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F)

plot(MainModel)

## Save image with results, currently commented to avoid accidental overwriting
#save.image("MainResults.RData")