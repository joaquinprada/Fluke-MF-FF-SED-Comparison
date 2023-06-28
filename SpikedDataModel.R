#########################
### Spiked Data model
### By JM Prada

## Load data and set working directory
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("SpikedInfection.csv")

## Change to real number of eggs counted
Flotac <- dt$Mini.FLOTAC...EPG/5 #0.2g sample
FKF <- dt$Flukefinder..EPG*2 #2g sample
Sed <- dt$Sedimentation...EPG*10 #10g sample

## model to estimate dose-sample relationship
nSamples <- nrow(dt)
Dosage <- dt$Dose*10 #Change EPGs to biggest sample - 10g
Species <- as.factor(dt$Parasite)

## Set seed ##
.RNG.seed <- function(chain)
  return( switch(chain, "1"= 1, "2"= 2) )
.RNG.name <- function(chain)
  return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )

## Tried with species specific rates, but both species give the same posterior
## No difference by species (so species removed from model)
m <- "model {
  
  for (n in 1:nSamples){
    
    ## Flotac
    Flotac[n] ~ dnegbin( (rt_Ft/(Dosage[n]/50+rt_Ft)), rt_Ft)
    
    ## Flukefinder
    FKF[n] ~ dnegbin( (rt_FF/(Dosage[n]/5+rt_FF)), rt_FF)
    
    ## Sedimentation
    Sed[n] ~ dnegbin( (rt_Sd/(Dosage[n]+rt_Sd)), rt_Sd)
  }
  
  rt_Ft ~ dgamma(0.001,0.001)
  rt_FF ~ dgamma(0.001,0.001)
  rt_Sd ~ dgamma(0.001,0.001)
  
  #inits# .RNG.seed, .RNG.name,
  #data# nSamples, Flotac, FKF, Sed, Dosage
  #monitor# rt_Ft, rt_FF, rt_Sd
}"

### Run model ###
NoiseSpiked <- run.jags(m, burnin=1000, sample=10000, thin=10, n.chains=2, jags.refresh = 1, method = 'parallel',
                        plots = F, silent.jags = F)

plot(NoiseSpiked)

rm(list=setdiff(ls(), "NoiseSpiked"))

## Save image with results, currently commented to avoid accidental overwriting
#save.image("SpikedResults.RData")