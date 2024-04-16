# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs/NewDataset")


# Load formatted data
load("dataset/formattedDataList.RData")

#load the predictions of the lengthAtAge
load("result/fittedModelLengthAtAge.RData")

#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
#nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

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
  
############################
# We assume that spawning probability is constant across the year ans
# assume a prior for the spawning probability
###########################
  spawnProb ~ dunif(0.0001, 0.3)
  
  #############################
  # Fecundicity
  #############################
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
      fecundicity[ind, age] <- exp(log(exp(lengthAtAgeThisYear[ind, age]))*2.21 - 6.15) * spawnProb
    }
  }
  
  ##############################################   
  # survival probability
  ###########################################
  #prior distributions
  
  for(i in 1:22){
    bsurv[i] ~dnorm(0, sd = 10)
  }
  
  sdbsurvAge ~ dunif(0.01, 10)
  sdbsurvLake ~ dunif(0.01, 10)
  #sdbsurvSpecies ~ dunif(0.01, 10)
  sdbsurvSex ~ dunif(0.01, 10)
  
  for(i in 1:maxAge){
    bsurvAge[i]~dnorm(0, sd = sdbsurvAge)
  }
  
  #lake effect
  for(i in 1:nLakes){
    bsurvLake[i] ~dnorm(0, sd = sdbsurvLake)
  }
  
  # for(i in 1:nSpecies){
  #   bsurvSpecies[i] ~ dnorm(0, sd = sdbsurvSpecies)
  # }
  
  for(i in 1:nSex){
    bsurvSex[i] ~ dnorm(0, sd = sdbsurvSex)
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
      logit(survivalProb[ind, age]) <- bsurv[1] + bsurvAge[age] + bsurvLake[lakeID[ind]] + bsurvSex[sexID[ind]] +  bsurvSpecies[speciesID[ind]]+ bsurv[2]* year[ind]  + bsurv[3]*forest[ind] + bsurv[4]*pastures[ind] + bsurv[5]*popnDensity[ind] + bsurv[6]*summerTemp[ind] + bsurv[7]*summerPrec[ind] + + bsurv[8]*summerSnowDepth[ind]
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

data = list(
  #sprawns = as.numeric(dataList$sprawningData$maturation),
  spawns = spawnsData,
  ageAtYear = as.numeric(dataList$sprawningData$ageAtYear),
  lengthAtAgeThisYear = lengthAtAgePosteriorSummary$data,
  year = as.numeric(scale(dataList$ageAtHarvestData$year)),
  y = as.matrix(dataList$ageAtHarvestData[ ,4:17]),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4:17])),
 forest = as.numeric(scale(dataList$ageAtHarvestData$`Coniferous forest`)),
  moorsAndHeathland = as.numeric(scale(dataList$ageAtHarvestData$`Moors and heathland`)),
  peatBogs = as.numeric(scale(dataList$ageAtHarvestData$`Peat bogs`)),
  waterbodies = as.numeric(scale(dataList$ageAtHarvestData$`Water bodies`)),
  broadleavedForest = as.numeric(scale(dataList$ageAtHarvestData$`Broad-leaved forest`)),
  sparselyVegAreas = as.numeric(scale(dataList$ageAtHarvestData$`Sparsely vegetated areas`)),
  transitionalWoodland = as.numeric(scale(dataList$ageAtHarvestData$`Transitional woodland-shrub`)),
  landOccByAgric = as.numeric(scale(dataList$ageAtHarvestData$`Land principally occupied by agriculture, with significant areas of natural vegetation`)),
  discontUrbanFabric = as.numeric(scale(dataList$ageAtHarvestData$`Discontinuous urban fabric`)),
  mixedForest = as.numeric(scale(dataList$ageAtHarvestData$`Mixed forest`)),
  seaAndOcean = as.numeric(scale(dataList$ageAtHarvestData$`Sea and ocean`)),
  bareRocks = as.numeric(scale(dataList$ageAtHarvestData$`Bare rocks`)),
  glaciersAndPerpetualSnow = as.numeric(scale(dataList$ageAtHarvestData$`Glaciers and perpetual snow`)),
  complexCultivationPatterns = as.numeric(scale(dataList$ageAtHarvestData$`Complex cultivation patterns`)),
  sportsAndLeisureFacilities = as.numeric(scale(dataList$ageAtHarvestData$`Sport and leisure facilities`)),
  greenUrbanAreas = as.numeric(scale(dataList$ageAtHarvestData$`Green urban areas`)),
  pastures = as.numeric(scale(dataList$ageAtHarvestData$Pastures)),
 summerTemp = as.numeric(scale(dataList$ageAtHarvestData$summerTemperature)),
 summerPrec = as.numeric(scale(dataList$ageAtHarvestData$summerPrecipitation)),
 summerSnowDepth = as.numeric(scale(dataList$ageAtHarvestData$summerSnowDepth)),
 popnDensity = as.numeric(scale(dataList$ageAtHarvestData$pop_density))
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sexID = as.numeric(dataList$sprawningData$sex),
  nSex = length(unique(as.numeric(dataList$sprawningData$sex))),
  lakeID = as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)),
  nLakes = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)))),
  maxAge = dim(data$y)[2],
  speciesID = as.numeric(as.factor(dataList$ageAtHarvestData$speciesID)),
  nSpecies = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$speciesID))))
)

inits <- list(
  bsp = rnorm(5, 0, 1),
  bl = rnorm(7, 0, 1),
  bsurv = rnorm(22, 0, 1),
  blSex = rnorm(constants$nSex, 0, 1),
  blSpecies = rnorm(constants$nSpecies, 0, 1),
  blAge = rnorm(constants$maxAge, 0, 1),
  blLake = rnorm(constants$nLakes, 0,1),
  bspSex = rnorm(constants$nSex, 0, 1),
  bspSpecies = rnorm(constants$nSpecies, 0, 1),
  bspAge = rnorm(constants$maxAge, 0, 1),
  bspLake = rnorm(constants$nLakes, 0,1),
  bsurvSex = rnorm(constants$nSex, 0, 1),
  bsurvSpecies = rnorm(constants$nSpecies, 0, 1),
  bsurvAge = rnorm(constants$maxAge, 0, 1),
  bsurvLake = rnorm(constants$nLakes, 0,1),
  sdblSex = 1,
  sdblSpecies = 1,
  sdblAge = 1,
  sdblLake = 1,
  sdbspSex = 1,
  sdbspSpecies = 1,
  sdbspAge = 1,
  sdbspLake = 1,
  sdbsurvSex = 1,
  sdbsurvSpecies = 1,
  sdbsurvAge = 1,
  sdbsurvLake = 1,
  Amat = array(0, dim = c(constants$maxAge,constants$maxAge, constants$nInds)),
  sdEpsilon = 1, 
  N = array(10, dim = c(constants$maxAge,constants$T, constants$nInds)),
  theta = matrix(0.1, nrow = constants$maxAge, ncol = constants$nInds),
  spawnProb = 0.1
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

# Configure the model
fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bsurv", "bsurvLake", "bsurvAge", "bsurvSpecies", "bsurvSex",
                                                  "sdbsurvLake", "sdbsurvAge",  "sdbsurvSex",
                                                  "lambda", "survivalProb","spawnProb")
)

# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
fishModelMCMCrun <- runMCMC(fishModelCompiled2,
                            niter = 300,
                            nchains = 2,
                            nburnin = 100,
                            thin = 2,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
fishModelMCMCrun$summary$all.chains

save(fishModelMCMCrun, file = "result/fittedModel.RData")

# currently the model runs. Will have to think about the model structure in details.