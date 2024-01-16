# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs")


# Load formatted data
load("dataset/formattedDataList.RData")


# Load required packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)


# start writing the nimble code to fit the model

code <- nimbleCode({
  #######################################
  # prior distributions
  ###########################################
  
  ## sprawing probability
  
  #covariate effect
 for(i in 1:4){
   bs[i] ~ dnorm(0, sd=100)
 }
  
  sdbs ~dunif(0.01, 100)
  
#treat sex as a factor in the sprawint probability  
  for(i in 1:2){
    bsSex[i] ~ dnorm(0, sd=sdbs)}
  
#variable selection Probability
  for(k in 1:5){
    psibs[k] ~ dunif(0,1)
  }

  # Sprawning probability
  for(ind in 1:nInds){
    logit(sprawnProb[ind]) <- bs[1] + psibs[1]*bs[2]*ageAtYear[ind] +  psibs[2]*bs[3]*lengthAtAgeThisYear[ind]  + psibs[3]*bs[3]*ageAtYear[ind]*lengthAtAgeThisYear[ind] + psibs[4]*bs[4]*CPUE[ind] +psibs[5]* bsSex[sex[ind]]
    }
 
  for(ind in 1:nInds){
    sprawns[ind] ~ dbern(sprawnProb[ind])
  }
  

   
}
)



# Data for model
data = list(
  sprawns = as.numeric(dataList$sprawningData$maturation),
  ageAtYear = as.numeric(dataList$sprawningData$ageAtYr),
  lengthAtAgeThisYear = as.numeric(dataList$sprawningData$legthAtAgeThisYear),
  CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort)
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sex = as.numeric(dataList$sprawningData$sex)
  )

inits <- list(
  bs = rnorm(4, 0, 1),
  bsSex = rnorm(2, 0, 1),
  sdbs = 1,
  psibs = rep(0.5, 5)
)



# nimbleModel
fishModel <- nimbleModel(code,
                         data = data,
                         constants = constants,
                         inits = inits)

# compile nimbleModel
fishModelCompiled <- compileNimble(fishModel)


# Configure fishModel
fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bs", "psibs", "sdbs")
)


# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
fishModelMCMCrun <- runMCMC(fishModelCompiled2,
                            niter = 1000,
                            nchains = 2,
                            nburnin = 100,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
fishModelMCMCrun$summary$all.chains
