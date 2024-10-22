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

# Write nimble Function to calculate eigen values and extract the maximum
lambdaEstimation <- function(x){
  ret <- max(Re(eigen(x)$values))
  return(ret)
}

nimbleLambdaEstimation <- nimble::nimbleRcall(
  prototype = function(x = double(2)){},
  returnType = double(0),
  Rfun = 'lambdaEstimation'
)

#nimbleLambdaEstimation(fishModelCompiled$Amat[,,1])


# start writing the nimble code to fit the model

code <- nimbleCode({
  
  #############################
  # Length at age this year
  #################################
  sdEpsilon ~ dunif(0.01, 10)
  
  for(i in 1:7){
    bl[i] ~ dnorm(0, sd = 10)
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){
  lengthLinearPred[ind, age] <- bl[1] + bl[2]*temp[ind] + bl[3]*precMay[ind] + bl[4]*precJune[ind] + bl[5]*precJuly[ind] + bl[6]*CPUE[ind] + bl[7]*NAO[ind]
  lengthAtAgeThisYear[ind, age] ~ dnorm(lengthLinearPred[ind, age], sd = sdEpsilon)
  }
  }
  
  #######################################
  # Sprawning probability
  ###########################################
  
  ## Prior distributions
  #covariate effect
  for(i in 1:5){
    bs[i] ~ dnorm(0, sd=10)
  }
  sdbs ~dunif(0.01, 10)
  sdbsAge ~dunif(0.01, 10)
  #treat sex as a factor in the sprawing probability  
  for(i in 1:2){
    bsSex[i] ~ dnorm(0, sd=sdbs)
  }
  
  #treat each stage/age category as a factor for the spawning probability
  for(age in 1:maxAge){
    bsAge[age] ~ dnorm(0, sd=sdbsAge)
  }
  
  # Add lake effect
  #lake effect
  sdLake ~ dunif(0.01, 10)
  for(l in 1:nLakes){
    bsLake[l] ~ dnorm(0, sd = sdLake)
  }
  
  #variable selection Probability
  for(k in 1:5){
    psibs[k] ~ dunif(0,1)
  }
  

  # Sprawning probability
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
    logit(sprawnProb[ind, age]) <- bs[1] + bsAge[age] + psibs[1]*bs[2]*ageAtYear[ind] +  psibs[2]*bs[3]*lengthAtAgeThisYear[ind]  + psibs[3]*bs[4]*ageAtYear[ind]*lengthAtAgeThisYear[ind] + psibs[4]*bs[5]*CPUE[ind] +psibs[5]* bsSex[sex[ind]] + bsLake[lake[ind]]
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
    fecundicity[ind, age] <- exp(log(lengthAtAgeThisYear[ind, age])*2.21 - 6.15) * sprawnProb[ind, age]
    }
  }
 
##############################################   
# survival probability
###########################################
 #prior distributions

  for(i in 1:9){
    bsurv[i] ~dnorm(0, sd = 10)
  }
  
  sdSurvAge ~ dunif(0.01, 10)
  sdSurvLake ~ dunif(0.01, 10)
  
  for(i in 1:maxAge){
    bsurvAge[i]~dnorm(0, sd = sdSurvAge)
  }
  
  #lake effect
  for(i in 1:nLakes){
   bsSurvLake[i] ~dnorm(0, sd = sdSurvLake)
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
 logit(survivalProb[ind, age]) <- bsurv[1] + bsurvAge[age] + bsurv[2]* yearGrowthOcc[ind] + bsSurvLake[lake[ind]] + bsurv[3]*forest[ind] + bsurv[4]*moorsAndHeathland[ind] + bsurv[5]*peatBogs[ind] + bsurv[6]*waterbodies[ind] + bsurv[7]*broadleavedForest[ind] + bsurv[8]*sparselyVegAreas[ind] + bsurv[9]*meanDVI[ind]  
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
    Amat[maxAge, maxAge, ind] <- survivalProb[ind, maxAge]
#  }

# Define population size
  #  for(ind in 1:nInds){
   #   N[1:maxAge, 1, ind] <- Nst[1:maxAge]
   #   for(t in 2:T){    
    #Population projection
  #  N[1:maxAge, t, ind] <- Amat[1:maxAge, 1:maxAge,ind]%*%N[1:maxAge, t-1, ind]
    
    #annual growth rate on log-scale
   # r.annual[t, ind] <- log(sum(N[1:maxAge, t, ind])) - log(sum(N[1:maxAge, t-1, ind]))
  }
  #} 
  
  # Growth rate (lambda) for each individual
  for(ind in 1:nInds){
   #lambda[ind] <- exp(mean(r.annual[u:T, ind])) 
    lambda[ind] <- nimbleLambdaEstimation(Amat[1:maxAge, 1:maxAge, ind])
    # Proportion of stable distribution
   # check skelly et. al (2023) page 1066 for equations
   
   # w = proportion to stable age
   # z = probability of remaining in the same age class
   # g = probability of surving to the next stage class
   w[ind, 1] <- 1 
   g[ind, 1] <- Amat[1,1,ind]
   g[ind, maxAge] <- 0
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
  for(age in 2:(maxAge)){
    ey[ind, age] <- g[ind, age-1]*C[ind, age-1] + z[ind, age]*C[ind, age]
  }
  #ey[ind, maxAge] <- g[ind, maxAge-1]*C[ind, maxAge-1] + survivalProb[ind, maxAge]*C[ind, maxAge]
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
  lengthAtAgeThisYear = as.matrix(dataList$lengthAtAgeThisYear[ ,4: 13]),#as.numeric(dataList$sprawningData$legthAtAgeThisYear),
  CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort),
  yearGrowthOcc = as.numeric(as.factor(dataList$ageAtHarvestData$yearGrowthOcc)),
  y = as.matrix(dataList$ageAtHarvestData[ ,4: 13]),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4: 13])),
  temp = as.numeric(scale(dataList$ageAtHarvestData$temperatureMeanInJuly)),
  precMay = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInMay)),
  precJune = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInJune)),
  precJuly = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInJuly)),
  forest = as.numeric(scale(dataList$ageAtHarvestData$`Coniferous forest`)),
  moorsAndHeathland = as.numeric(scale(dataList$ageAtHarvestData$`Moors and heathland`)),
  peatBogs = as.numeric(scale(dataList$ageAtHarvestData$`Peat bogs`)),
  waterbodies = as.numeric(scale(dataList$ageAtHarvestData$`Water bodies`)),
  broadleavedForest = as.numeric(scale(dataList$ageAtHarvestData$`Broad-leaved forest`)),
  sparselyVegAreas = as.numeric(scale(dataList$ageAtHarvestData$`Sparsely vegetated areas`)),
  meanDVI = as.numeric(dataList$ageAtHarvestData$mean_ndvi),
  NAO = as.numeric(dataList$ageAtHarvestData$NAOwinterIndex)
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sex = as.numeric(dataList$sprawningData$sex),
  lake = as.numeric(as.factor(dataList$ageAtHarvestData$lakeName)),
  nLakes = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$lakeName)))),
  maxAge = 10,
  Nst = rep(3, 10),
  T = 100,
  u = 80 # cut off year from the prediction of population size used to estimate popn growth rate
)

inits <- list(
  bs = rnorm(5, 0, 1),
  bl = rnorm(7, 0, 1),
  bsSex = rnorm(2, 0, 1),
  bsAge = rnorm(constants$maxAge, 0, 1),
  sdbsAge =1,
  sdLake = 1,
  sdSurvLake = 1,
  bsLake = rnorm(constants$nLakes, 0, 1),
  bsSurvLake = rnorm(constants$nLakes, 0, 1),
  sdbs = 1,
  psibs = rep(1, 5),
  bsurv = rnorm(9, 0, 1),
  bsurvAge = rnorm(constants$maxAge, 0, 1),
  sdSurvAge = 1,
  Amat = array(0, dim = c(10,10, constants$nInds)),
  sdEpsilon = 1, 
  N = array(10, dim = c(10,constants$T, constants$nInds)),
  theta = matrix(0.1, nrow = constants$maxAge, ncol = constants$nInds)
)



# nimbleModel
fishModel <- nimbleModel(code,
                         data = data,
                         constants = constants,
                         inits = inits)

# compile nimbleModel
fishModelCompiled <- compileNimble(fishModel)

# check if any lambda is 0
which(fishModelCompiled$lambda == 0)

# Configure fishModel
fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bs", "psibs", "sdbs",
                                                  "bsSex", "bsSurvLake", "bsLake",
                                                  "bsurvAge", "bsurv",
                                                  "sdbsAge", "bsAge",
                                                   "sdEpsilon", "lambda", "bl", "survivalProb")
)
#"lambda",

# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
fishModelMCMCrun <- runMCMC(fishModelCompiled2,
                            niter = 2000,
                            nchains = 2,
                            nburnin = 1000,
                            thin = 2,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
fishModelMCMCrun$summary$all.chains

save(fishModelMCMCrun, file = "modelFitting/fittedModel.RData")

# currently the model runs. Will have to think about the model structure in details.