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
#  for(i in 1:4){
#    bs[i] ~ dnorm(0, sd=100)
#  }
#   
#   sdbs ~dunif(0.01, 100)
#   
# #treat sex as a factor in the sprawint probability  
#   for(i in 1:2){
#     bsSex[i] ~ dnorm(0, sd=sdbs)}
#   
# #variable selection Probability
#   for(k in 1:5){
#     psibs[k] ~ dunif(0,1)
#   }
# 
#   # Sprawning probability
#   for(ind in 1:nInds){
#     logit(sprawnProb[ind]) <- bs[1] + psibs[1]*bs[2]*ageAtYear[ind] +  psibs[2]*bs[3]*lengthAtAgeThisYear[ind]  + psibs[3]*bs[3]*ageAtYear[ind]*lengthAtAgeThisYear[ind] + psibs[4]*bs[4]*CPUE[ind] +psibs[5]* bsSex[sex[ind]]
#     }
#  
#   for(ind in 1:nInds){
#     sprawns[ind] ~ dbern(sprawnProb[ind])
#   }
#   

  for(i in 1:4){
    bs[i] ~ dnorm(0, sd=100)
  }
  sdbs ~dunif(0.01, 100)
  sdbsAge ~dunif(0.01, 100)
  #treat sex as a factor in the sprawing probability  
  for(i in 1:2){
    bsSex[i] ~ dnorm(0, sd=sdbs)
  }
  
  #treat each stage/age category as a factor for the spawning probability
  for(age in 1:maxAge){
    bsAge[age] ~ dnorm(0, sd=sdbsAge)
  }
  
  #variable selection Probability
  for(k in 1:5){
    psibs[k] ~ dunif(0,1)
  }
  
  
  # Sprawning probability
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
      logit(sprawnProb[ind, age]) <- bs[1] + bsAge[age] + psibs[1]*bs[2]*ageAtYear[ind] +  psibs[2]*bs[3]*lengthAtAgeThisYear[ind]  + psibs[3]*bs[3]*ageAtYear[ind]*lengthAtAgeThisYear[ind] + psibs[4]*bs[4]*CPUE[ind] +psibs[5]* bsSex[sex[ind]]
    }
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
      sprawns[ind, age] ~ dbern(sprawnProb[ind, age])
    }
  }
   
}
)



# # Data for model
# data = list(
#   sprawns = as.numeric(dataList$sprawningData$maturation),
#   ageAtYear = as.numeric(dataList$sprawningData$ageAtYr),
#   lengthAtAgeThisYear = as.numeric(dataList$sprawningData$legthAtAgeThisYear),
#   CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort)
# )
# 
# constants = list(
#   nInds = nrow(dataList$ageAtHarvestData),
#   sex = as.numeric(dataList$sprawningData$sex)
#   )
# 
# inits <- list(
#   bs = rnorm(4, 0, 1),
#   bsSex = rnorm(2, 0, 1),
#   sdbs = 1,
#   psibs = rep(0.5, 5)
# )
# Data for model

sprawnsData <- apply(dataList$sprawningData[, 4:13], c(1,2), function(x) ifelse(x>0, 1, 0))

data = list(
  #sprawns = as.numeric(dataList$sprawningData$maturation),
  sprawns = sprawnsData,
  ageAtYear = as.numeric(dataList$sprawningData$ageAtYr),
  lengthAtAgeThisYear = as.numeric(dataList$sprawningData$legthAtAgeThisYear),
  CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort),
  yearGrowthOcc = as.numeric(as.factor(dataList$ageAtHarvestData$yearGrowthOcc)),
  y = as.matrix(dataList$ageAtHarvestData[ ,4: 13]),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4: 13]))
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sex = as.numeric(dataList$sprawningData$sex),
  maxAge = 10,
  Nst = rep(3, 10),
  T = 100,
  u = 80 # cut off year from the prediction of population size used to estimate popn growth rate
)

inits <- list(
  bs = rnorm(4, 0, 1),
  bsSex = rnorm(2, 0, 1),
  bsAge = rnorm(constants$maxAge, 0, 1),
  sdbsAge =1,
  sdbs = 1,
  psibs = rep(0.5, 5),
  bsurv = rnorm(2, 0, 1),
  bsurvAge = rnorm(constants$maxAge, 0, 1),
  sdSurvAge = 1,
  Amat = array(0, dim = c(10,10, constants$nInds)),
  N = array(10, dim = c(10,constants$T, constants$nInds)),
  theta <- matrix(0.1, nrow = constants$maxAge, ncol = constants$nInds)
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
                                     monitors = c("bs", "psibs", "sdbs",
                                                  "bsSex", "bsAge")
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
