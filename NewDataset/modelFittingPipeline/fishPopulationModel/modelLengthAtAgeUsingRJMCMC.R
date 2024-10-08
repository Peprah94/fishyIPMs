# This script is used to model the length at age of the fish species
# and predict the lengths that are missing
# We also use the RJMCMC to perform variable selection of the covariates
# Read the RJMCMC documentation here: https://r-nimble.org/writing-reversible-jump-mcmc-in-nimble


# set working directory.
#change this for your particular application
#setwd("C:/GitHub/fishyIPMs/NewDataset")

set.seed(1994)

# Load formatted data
load("dataset/formattedDataList.RData")


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
  #sdEpsilon ~ dunif(0.01, 10)
 # sdblLake ~ dunif(0.01, 10)
 # sdblAge ~ dunif(0.01, 10)
  
  
  #sdblSex ~ dunif(0.01, 10)
 
  sdEpsilon ~ dunif(0.01, 10)
   sdblLake <- 10
   sdblAge <- 10
  sdblSex <- 10
  
  intercept ~ dnorm(0, sd = 10)
  
  
  # for(i in 1:nLakes){
  #   blLake[i] ~ dnorm(0, sd = sdblLake)
  # }
  # 
  # for(i in 1:nSex){
  #   blSex[i] ~ dnorm(0, sd = sdblSex)
  # }
  # 
   for(i in 1:maxAge){
     blAge[i] ~ dnorm(0, sd= sdblAge)
   }
  # 
  # Inclusion of variables in a model
  alpha ~ dunif(0,1)    ## prior on inclusion probability

  
  for(i in 1:12){
  pbLength[i] ~ dbern(alpha)  ## indicator variable for each coefficient
  bl[i] ~ dnorm(0, 10)
  blIndex[i] <- bl[i] * pbLength[i]
  }
  
  # + blSpecies[speciesID[ind]]
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){
      lengthLinearPred[ind, age] <- blAge[age] + intercept + 
                                    blIndex[1]*year[ind] + 
                                    blIndex[2]*summerTemp[ind] + 
                                    blIndex[3]*summerSnowDepth[ind] + 
                                    blIndex[4]*forest[ind] + 
                                    blIndex[5]*WPUE[ind] + 
                                    blIndex[6]*sexID[ind]  + 
                                    blIndex[7]*lakeID[ind] +
                                    blIndex[8]* pow(summerTemp[ind], 2) +
                                    blIndex[9]* pow(summerSnowDepth[ind], 2)+
                                    blIndex[10] * summerTemp[ind] * summerSnowDepth[ind]+
                                    blIndex[11] * summerTemp[ind] * WPUE[ind] + 
                                    blIndex[12] * summerSnowDepth[ind] * WPUE[ind] 
      #lengthLinearPred[ind, age] <- blAge[age]  + blSex[sexID[ind]] + blLake[lakeID[ind]] + bl[1] + pbLength[1]*bl[2]*year[ind] + pbLength[2]*bl[3]*summerTemp[ind] + pbLength[3]*bl[4]*summerSnowDepth[ind] +pbLength[4]* bl[5]*forest[ind]
      
      lengthAtAgeThisYear[ind, age] ~ dnorm(lengthLinearPred[ind, age], sd = sdEpsilon)
    }
  }
  
  
  
}
)



# Data for model
stagedData <- rowSums(as.matrix(dataList$ageAtHarvestData[ ,9:17])) %>%
  cbind(as.matrix(dataList$ageAtHarvestData[ ,4:8]), .)

stagedLengthAtAge <- rowMeans(as.matrix(dataList$lengthAtAgeData[ ,9: 17]), 
                              na.rm = TRUE)%>%
                      cbind(as.matrix(dataList$lengthAtAgeData[ ,4: 8]), .)
stagedLengthAtAge[is.nan(stagedLengthAtAge)] <- NA

#spawnsData <- apply(dataList$sprawningData[, 4:17], c(1,2), function(x) ifelse(x>0, 1, 0))
#wpue = as.numeric(scale(dataList$ageAtHarvestData$WPUE))
#wpue[is.na(wpue)] <- mean(wpue, na.rm = TRUE)

data = list(
  #sprawns = as.numeric(dataList$sprawningData$maturation),
  #spawns = spawnsData,
  ageAtYear = as.numeric(dataList$sprawningData$ageAtYear),
  lengthAtAgeThisYear = log(stagedLengthAtAge),#as.numeric(dataList$sprawningData$legthAtAgeThisYear),
  # CPUE = as.numeric(dataList$sprawningData$capturePerUnitEffort),
  year = as.numeric(scale(dataList$ageAtHarvestData$year)),#as.numeric(as.factor(dataList$ageAtHarvestData$yearGrowthOcc)),
  y = as.matrix(stagedData),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4:17])),
  lakeID = (as.numeric(scale(dataList$ageAtHarvestData$Y))), 
  forest = as.numeric(scale(dataList$ageAtHarvestData$forest)),
  summerTemp = as.numeric(scale(dataList$ageAtHarvestData$meanSummerTemp)),
  summerSnowDepth = as.numeric(scale(dataList$ageAtHarvestData$meanWinterSnow)),
  popnDensity = as.numeric(scale(dataList$ageAtHarvestData$popnDensity)),
  WPUE = as.numeric(scale(dataList$ageAtHarvestData$WPUE))
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sexID = (as.numeric(dataList$sprawningData$sex)-1), # Make one category the reference category
  nSex = length(unique(as.numeric(dataList$sprawningData$sex))),
  nLakes = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)))),
  maxAge = dim(data$y)[2],
  speciesID = as.numeric(as.factor(dataList$ageAtHarvestData$speciesID)),
  nSpecies = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$speciesID))))
)

inits <- list(
  pbLength = rep(1, 12),
  bsp = rnorm(5, 0, 1),
  bl = rnorm(12, 0, 1),
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

# Configure model

fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bl", "blAge", 
                                                  "sdblLake", "sdblAge", "sdblSex", 
                                                  "lengthAtAgeThisYear", "sdEpsilon",
                                                  "pbLength", "alpha", "intercept")
)

#"blLake", "blSex",
configureRJ(fishModelConfigured,
            targetNodes = c("bl"),
            indicatorNodes = 'pbLength',
            control = list(mean = 0, scale = .2))

fishModelConfigured$printSamplers(
)
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

save(lengthAtAgePosteriorSummary, 
     file = "result/fittedModelLengthAtAge.RData")

save(lengthAtAgeModelRun,
     file = "result/fittedModelLengthAtAgeMCMCRun.RData")
# currently the model runs. Will have to think about the model structure in details.