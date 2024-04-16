# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs/NewDataset")


# Load formatted data
load("dataset/formattedDataList.RData")

#nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
#nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE)

# Load required packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)


# start writing the nimble code to fit the model

code <- nimbleCode({
  
  #############################
  # Length at age this year
  #################################
  sdEpsilon ~ dunif(0.01, 10)
  sdblLake ~ dunif(0.01, 10)
  sdblAge ~ dunif(0.01, 10)
  #sdblSpecies ~ dunif(0.01, 10)
  sdblSex ~ dunif(0.01, 10)
  for(i in 1:6){
    bl[i] ~ dnorm(0, sd = 10)
  }
  
  for(i in 1:nLakes){
    blLake[i] ~ dnorm(0, sd = sdblLake)
  }
  
  for(i in 1:nSex){
    blSex[i] ~ dnorm(0, sd = sdblSex)
  }
  
  for(i in 1:maxAge){
    blAge[i] ~dnorm(0, sd= sdblAge)
  }
  
  # for(i in 1:nSpecies){
  #   blSpecies[i] ~ dnorm(0, sdblSpecies)
  # }
  
  # + blSpecies[speciesID[ind]]
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){
      lengthLinearPred[ind, age] <- blAge[age]  + blSex[sexID[ind]] + blLake[lakeID[ind]] + bl[1] + bl[2]*year[ind] + bl[3]*summerTemp[ind] + bl[4]*summerPrec[ind] + bl[5]*summerSnowDepth[ind] + bl[6]*forest[ind]
      #lengthLinearPred[ind, age] <- bl[1] + bl[2]*temp[ind] + bl[3]*precMay[ind] + bl[4]*precJune[ind] + bl[5]*precJuly[ind] + bl[6]*CPUE[ind] + bl[7]*NAO[ind]
      lengthAtAgeThisYear[ind, age] ~ dnorm(lengthLinearPred[ind, age], sd = sdEpsilon)
    }
  }
  
  
  
}
)



# Data for model

spawnsData <- apply(dataList$sprawningData[, 4:17], c(1,2), function(x) ifelse(x>0, 1, 0))


data = list(
  #sprawns = as.numeric(dataList$sprawningData$maturation),
  spawns = spawnsData,
  ageAtYear = as.numeric(dataList$sprawningData$ageAtYear),
  lengthAtAgeThisYear = log(as.matrix(dataList$lengthAtAgeData[ ,4: 17])),#as.numeric(dataList$sprawningData$legthAtAgeThisYear),
  # CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort),
  year = as.numeric(scale(dataList$ageAtHarvestData$year)),#as.numeric(as.factor(dataList$ageAtHarvestData$yearGrowthOcc)),
  y = as.matrix(dataList$ageAtHarvestData[ ,4:17]),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4:17])),
  #temp = as.numeric(scale(dataList$ageAtHarvestData$temperatureMeanInJuly)),
  # precMay = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInMay)),
  #precJune = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInJune)),
  # precJuly = as.numeric(scale(dataList$ageAtHarvestData$precipitationSumInJuly)),
  forest = as.numeric(scale(dataList$ageAtHarvestData$`Coniferous forest`)),
  moorsAndHeathland = as.numeric(scale(dataList$ageAtHarvestData$`Moors and heathland`)),
  peatBogs = as.numeric(scale(dataList$ageAtHarvestData$`Peat bogs`)),
  waterbodies = as.numeric(scale(dataList$ageAtHarvestData$`Water bodies`)),
  broadleavedForest = as.numeric(scale(dataList$ageAtHarvestData$`Broad-leaved forest`)),
  sparselyVegAreas = as.numeric(scale(dataList$ageAtHarvestData$`Sparsely vegetated areas`)),
  #meanDVI = as.numeric(dataList$ageAtHarvestData$mean_ndvi),
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
  # NAO = as.numeric(dataList$ageAtHarvestData$NAOwinterIndex)
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
  bl = rnorm(6, 0, 1),
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
  sprawnProb = 0.1
)



# nimbleModel
fishModel <- nimbleModel(code,
                         data = data,
                         constants = constants,
                         inits = inits)

# compile nimbleModel
fishModelCompiled <- compileNimble(fishModel)

# Configure fishModel
# fishModelConfigured <- configureMCMC(fishModelCompiled,
#                                      monitors = c("bl", "blLake", "blAge", "blSpecies", "blSex",
#                                                   "bsp", "bspLake", "bspAge", "bspSpecies", "bspSex",
#                                                   "bsurv", "bsurvLake", "bsurvAge", "bsurvSpecies", "bsurvSex",
#                                                    "sdblLake", "sdblAge", "sdblSpecies", "sdblSex",
#                                                    "sdbspLake", "sdbspAge", "sdbspSpecies", "sdbspSex",
#                                                    "sdbsurvLake", "sdbsurvAge", "sdbsurvSpecies", "sdbsurvSex",
#                                                   "sdEpsilon", "lambda", "survivalProb","spawnProb")
# )

fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bl", "blLake", "blAge", "blSex",
                                                  "sdblLake", "sdblAge", "sdblSex", "lengthAtAgeThisYear")
)
#"lambda",

#change samplers for some MCMC configurations
#fishModelConfigured$removeSamplers(c("blLake", "blAge"))
#fishModelConfigured$addSampler(target = "blLake", type = "RW_block")
#fishModelConfigured$addSampler(target = "blAge", type = "RW_block")

# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
lengthAtAgeModelRun <- runMCMC(fishModelCompiled2,
                            niter = 30000,
                            nchains = 2,
                            nburnin = 10000,
                            thin = 2,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
lengthAtAgePosteriorSummary <- list()
lengthAtAgePosteriorSummary$data <- matrix(lengthAtAgeModelRun$summary$all.chains[rownames(lengthAtAgeModelRun$summary$all.chains)[grepl("lengthAtAge", rownames(lengthAtAgeModelRun$summary$all.chains))] , "Mean"],
                                           nrow = constants$nInds,
                                           ncol = constants$maxAge, 
                                           byrow = FALSE)
lengthAtAgePosteriorSummary$summary <- lengthAtAgeModelRun$summary$all.chains[rownames(lengthAtAgeModelRun$summary$all.chains)[!grepl("lengthAtAge", rownames(lengthAtAgeModelRun$summary$all.chains))] , ]

save(lengthAtAgePosteriorSummary, file = "result/fittedModelLengthAtAge.RData")

# currently the model runs. Will have to think about the model structure in details.