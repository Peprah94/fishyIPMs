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
  # Sprawning probability
  ###########################################
  
  ## Prior distributions
  #covariate effect
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
  
# Fecundicity
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
    fecundicity[ind, age] <- exp(log(lengthAtAgeThisYear[ind])*2.21 - 6.15) * sprawnProb[ind, age]
    }
  }
 
##############################################   
# survival probability
###########################################
 #prior distributions

  for(i in 1:2){
    bsurv[i] ~dnorm(0, sd = 100)
  }
  
  sdSurvAge ~dunif(0.01, 100)
  
  for(i in 1:maxAge){
    bsurvAge[i]~dnorm(0, sd = sdSurvAge)
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
 logit(survivalProb[ind, age]) <- bsurv[1] + bsurvAge[age] + bsurv[2]* yearGrowthOcc[ind]     
    }
  }
  
  ## Population size
  # Initialise the population size nodes
  #Amat[1:10, 1:10 ,1:nInds] <- 0
  for(ind in 1:nInds){
    Amat[1,1 ,ind] <-  fecundicity[ind, 1]*survivalProb[ind, 1]
    for(age in 2:maxAge){ #loop through the age
    Amat[1,age ,ind] <-  fecundicity[ind, age]*survivalProb[ind, 1]
    Amat[age,(age-1),ind] <-  survivalProb[ind, age]
    }
    Amat[maxAge, maxAge, ind] <- survivalProb[maxAge, maxAge]
#  }

# Define population size

  #  for(ind in 1:nInds){
      N[1:maxAge, 1, ind] <- Nst[1:maxAge]
      for(t in 2:T){    
    #Population projection
    N[1:maxAge, t, ind] <- Amat[1:maxAge, 1:maxAge,ind]%*%N[1:maxAge, t-1, ind]
    
    #annual growth rate on log-scale
    r.annual[t, ind] <- log(sum(N[1:maxAge, t, ind])) - log(sum(N[1:maxAge, t-1, ind]))
  }
  } 
  
  # Growth rate (lambda) for each individual
  for(ind in 1:nInds){
   lambda[ind] <- exp(mean(r.annual[u:T, ind])) 
  
    # Proportion of stable distribution
   # check skelly et. al (2023) page 1066 for equations
   
   # w = proportion to stable age
   # z = probability of remaining in the same age class
   # g = probability of surving to the next stage class
   w[ind, 1] <- 1 
   g[ind, 1] <- Amat[1,1,ind]
   g[ind, maxAge] <- 0.0000001
   z[ind, 1] <- 0
   z[ind, maxAge] <- Amat[maxAge,maxAge - 1,ind]
   for(age in 2:(maxAge-1)){
     z[ind, age] <- 0
     g[ind, age] <- Amat[age,(age-1),ind]
     w[ind, age] <- (g[ind, age-1]/(lambda[ind] -z[ind, age]))*w[ind, age-1]
   }
   w[ind, maxAge] <- (g[ind, maxAge]/(lambda[ind] -z[ind, maxAge]))*w[ind, maxAge-1] 
   
   #stable state distribution
  C[ind, 1:maxAge] <- w[ind, 1:maxAge]/sum(w[ind, 1:maxAge]) 
   
 
  #expected counts in each stage class
  for(age in 2:(maxAge-1)){
    ey[ind, age] <- g[ind, age-1]*C[ind, age-1] + z[ind, age]*C[ind, age]
  }
  ey[ind, maxAge] <- g[ind, maxAge-1]*C[ind, maxAge-1] + survivalProb[ind, maxAge]*C[ind, maxAge]
  ey[ind, 1] <- (lambda[ind] - sum(ey[ind, 2:maxAge]))
  
  
  # likelihood
  #theta[ind, 1:maxAge] ~ ddirch(ey[ind, 1:maxAge])
  y[ind, 1:maxAge] ~ dmulti(ey[ind, 1:maxAge], harvestCount[ind])
  
  }

}
)



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
                                                  "bsSex",
                                                  "bsurvAge", "bsurv",
                                                  "sdbsAge", "bsAge")
)


# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
fishModelMCMCrun <- runMCMC(fishModelCompiled2,
                            niter = 50000,
                            nchains = 2,
                            nburnin = 10000,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
fishModelMCMCrun$summary$all.chains

# currently the model runs. Will have to think about the model structure in details.