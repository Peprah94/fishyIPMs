# Simulation study

# set working directory.
#change this for your particular application
#setwd("C:/GitHub/fishyIPMs/NewDataset")
# Load required packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)

constrainedLambda <- function(z, lambda){
  # Get the maximum of z
  zMax <- max(z[1:length(z)])
  
  if(zMax >= lambda){
    lambdaInd <- zMax + 0.0001
  } else {
    lambdaInd <- lambda
  } 
  return(lambdaInd)
}

nimbleConstrainedLambda <- nimble::nimbleRcall(
  prototype = function(z = double(1), lambda = double(0)){},
  returnType = double(0),
  Rfun = 'constrainedLambda'
)

# Load formatted data
load("dataset/formattedDataList.RData")

#load the predictions of the lengthAtAge
load("result/fittedModelLengthAtAge.RData")

dataList$ageAtHarvestDataObs <- dataList$ageAtHarvestData %>%
  arrange(year)%>%
  filter(year < 2013)

dataList$lengthAtAgeDataObs <- dataList$lengthAtAgeData %>%
  arrange(year)%>%
  filter(year < 2013)

dataList$ageAtHarvestDataPred <- dataList$ageAtHarvestData %>%
  arrange(year)%>%
  filter(year >= 2013)

dataList$lengthAtAgeDataPred <- dataList$lengthAtAgeData %>%
  arrange(year)%>%
  filter(year >= 2013)



# Data for model
stagedData <- rowSums(as.matrix(dataList$ageAtHarvestDataObs[ ,9:17])) %>%
  cbind(as.matrix(dataList$ageAtHarvestDataObs[ ,4:8]), .)

stagedLengthAtAge <- rowMeans(as.matrix(dataList$lengthAtAgeDataObs[ ,9: 17]), 
                              na.rm = TRUE)%>%
  cbind(as.matrix(dataList$lengthAtAgeDataObs[ ,4: 8]), .)
stagedLengthAtAge[is.nan(stagedLengthAtAge)] <- NA

# Extract the covariates needed
forest <- as.numeric(scale(dataList$ageAtHarvestData$forest, scale = TRUE))
temp <- as.numeric(scale(dataList$ageAtHarvestData$meanSummerTemp, scale = TRUE))
year <- as.numeric(scale(dataList$ageAtHarvestData$year, scale = TRUE))
winterSnowDepth <- as.numeric(scale(dataList$ageAtHarvestData$meanWinterSnow, scale = TRUE))
popnDensity <- as.numeric(scale(dataList$ageAtHarvestData$popnDensity, scale = TRUE))
wpue <- as.numeric(scale(dataList$ageAtHarvestData$WPUE, scale = TRUE))
latitude <- as.numeric(scale(dataList$ageAtHarvestData$Y, scale = TRUE))

# Number od Obs and Predictions
nObs <- nrow(dataList$ageAtHarvestDataObs)
nPreds <- nrow(dataList$ageAtHarvestDataPred)

medLength <- apply(lengthAtAgePosteriorSummary$data, 2, median)

data = list(
  # ageAtYear = as.numeric(dataList$sprawningData$ageAtYear),
  lengthAtAgeThisYear = lengthAtAgePosteriorSummary$data,
  medLength = medLength,
  #year = as.numeric(as.factor(dataList$ageAtHarvestData$year)),
  y = as.matrix(stagedData),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestDataObs[4:17])),
  forest = forest[1:nObs],
  #pastures = as.numeric(scale(dataList$ageAtHarvestData$Pastures)),
  summerTemp = temp[1:nObs],
  lakeID = latitude[1:nObs], 
  #summerPrec = as.numeric(scale(dataList$ageAtHarvestData$summerPrecipitation)),
  summerSnowDepth = winterSnowDepth[1:nObs],
  popnDensity = popnDensity[1:nObs],
  forestsq = (forest[1:nObs])^2,
  # pasturessq = as.numeric((scale(dataList$ageAtHarvestData$Pastures)))^2,
  summerTempsq = (temp[1:nObs])^2,
  #summerPrecsq = as.numeric((scale(dataList$ageAtHarvestData$summerPrecipitation)))^2,
  summerSnowDepthsq = (winterSnowDepth[1:nObs])^2,
  popnDensitysq = (popnDensity[1:nObs])^2,
  WPUE = wpue[1:nObs],
  year = year[1:nObs],
  forestPred = forest[-c(1:nObs)],
  #pastures = as.numeric(scale(dataList$ageAtHarvestData$Pastures)),
  summerTempPred = temp[-c(1:nObs)],
  lakeIDPred = latitude[-c(1:nObs)], 
  #summerPrec = as.numeric(scale(dataList$ageAtHarvestData$summerPrecipitation)),
  summerSnowDepthPred = winterSnowDepth[-c(1:nObs)],
  popnDensityPred = popnDensity[-c(1:nObs)],
  forestsqPred = (forest[-c(1:nObs)])^2,
  # pasturessq = as.numeric((scale(dataList$ageAtHarvestData$Pastures)))^2,
  summerTempsqPred = (temp[-c(1:nObs)])^2,
  #summerPrecsq = as.numeric((scale(dataList$ageAtHarvestData$summerPrecipitation)))^2,
  summerSnowDepthsqPred = (winterSnowDepth[-c(1:nObs)])^2,
  popnDensitysqPred = (popnDensity[-c(1:nObs)])^2,
  WPUEPred = wpue[-c(1:nObs)],
  yearPred = year[-c(1:nObs)],
  d = c(rep(1, dim(as.matrix(stagedData))[2] - 1 ), 1)
)

constants = list(
  nInds = nObs,
  nPreds = nPreds, 
  sexIDObs = (as.numeric(dataList$ageAtHarvestData$sex))[1:nObs],
  sexIDPred = (as.numeric(dataList$ageAtHarvestData$sex))[-c(1:nObs)],
  nSex = length(unique(as.numeric(dataList$ageAtHarvestData$sex))),
  nLakes = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)))),
  maxAge = dim(data$y)[2]#,
  #,
  # nYears = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$year))))
)

inits <- list(
  psi = rep(1, 10),
  bsurv = rnorm(10, 0, 1),
  #bsurvIndex = rnorm(10, 0, 1),
  alpha = 0.5,
  intercept = 1,
  bsurvSex = rnorm(constants$nSex, 0, 1),
  bsurvAge = rnorm(constants$maxAge, 0, 1),
  bsurvProbSex = rnorm(constants$nSex, 0, 1),
  bsurvProbAge = runif(constants$maxAge, 0, 1),
  #bsurvYear = rnorm(constants$nYears, 0,1),
  sdbsurvSex = 1,
  sdbsurvAge = 1,
  sdbsurvYear = 1,
  sdbsurvProbSex = 1,
  sdbsurvProbAge = 1,
  sdbsurvProbYear = 1,
  survProbIntercept = 0, 
  bsurvProbYear = 0,
  #lambda = rep(1, constants$nYears),
  Amat = array(0, dim = c(constants$maxAge,constants$maxAge, constants$nInds)),
  spawnProb = 0.1,
  theta = matrix((rep(1, constants$maxAge)/constants$maxAge), nrow = constants$nInds, ncol = constants$maxAge)
)


simData <- function(data, constants, 
                    seed,
                    beta = c(- 1.5, - 0.5, 3, 0, 0, 0.8, -1), 
                    sdbsurvSex = 1,
                    sdbsurvAge = 1,
                    alpha = c(-1, -1.5, -0.5),
                    selectAlpha = 0.8){
  set.seed(seed)
  
  maxAge <- constants$maxAge
  d <- data$d
  # Simulate lambda
  lambda <- lambdaInd <-  vector("numeric", constants$nInds)
  bsurvSex <- vector("numeric", constants$nSex)
  survivalProb <- w <- g <- z <- y <- ey <- C <- theta <-selectivity <- matrix(NA, nrow = constants$nInds, ncol = constants$maxAge)
  
  # Simulate random effect for sex
    bsurvSex <- rnorm(constants$nSex,0, sdbsurvSex)
  
  
  # simulate random effect for age
    bsurvProbAge <- c(0.9, 0.5, 0.2, 0.1, 0.02, 0.01) #rnorm(constants$maxAge, 0, sdbsurvAge)  

  
  for(i in 1:constants$nInds){
  lambda[i] <- exp(beta[1] + 
                     beta[2]*data$year[i] + 
                     beta[3]*data$lakeID[i] + 
                     beta[4]* (data$year[i])^2 +
                     beta[5]* (data$lakeID[i])^2 +
                     beta[6]* (data$summerTemp[i]) +
                     beta[7] * (data$summerTemp[i])^2 +
                     bsurvSex[constants$sexIDObs[i]])
  }
  
  for(ind in 1:constants$nInds){
    for(age in 1:constants$maxAge){ 
      survivalProb[ind, age] <- bsurvProbAge[age]
                                        # alpha[2]*data$year[ind] +
                                         #alpha[3]* data$summerTemp[ind]+
                                         #bsurvProbYear *year[ind] + bsurvProbSex[sexIDObs[ind]] + bsurvProbAge[age] + survProbIntercept
    
      # Per capita fecundicity
      selectivity[ind, age] <- 1/(1+ exp(-selectAlpha * (data$lengthAtAgeThisYear[ind, age] - data$medLength[age])))
        
        #(exp(log(exp(data$lengthAtAgeThisYear[ind, age]))*2.21 - 6.15) * spawnProb)/data$harvestCount[ind]
      
      }
  }

  
    for(ind in 1:constants$nInds){
      #lambda[ind] <- exp(mean(r.annual[u:T, ind])) 
      
      # Asymptotic population growth rate
      
      
      #nimbleLambdaEstimation(Amat[1:maxAge, 1:maxAge, ind])
      # Proportion of stable distribution
      # check skelly et. al (2023) page 1066 for equations
      
      # w = proportion to stable age
      # z = probability of remaining in the same age class
      # g = probability of surviving to the next stage class
      w[ind, 1] <- 1 
      # Probability of surviving to the next stage class
      g[ind, 1] <- survivalProb[ind, 1]
      g[ind, maxAge] <- (survivalProb[ind, maxAge])^d[maxAge] * (1 - survivalProb[ind, maxAge])/ (1 - (survivalProb[ind, maxAge])^d[maxAge] )
      z[ind, maxAge] <- ((1 - survivalProb[ind, maxAge] ^ (d[maxAge]  - 1)))* survivalProb[ind, maxAge] / (1 - survivalProb[ind, maxAge] ^ d[maxAge] ) 
      z[ind, 1] <- 0 #survivalProb[ind, maxAge]
      for(age in 2:(maxAge-1)){
        z[ind, age] <- 0 #(1 - survivalProb[ind, age] ^ (d[age] - 1)) / (1 - survivalProb[ind, age] ^ d[age]) * survivalProb[ind, age]
        g[ind, age] <- (survivalProb[ind, age])^d[age] * (1 - survivalProb[ind, age])/ (1 - (survivalProb[ind, age])^d[age])
      }
      
      # Contrain the lambda estimates to avoid identifiability issues
      lambdaInd[ind] <- constrainedLambda(z[ind, 1:maxAge], lambda[ind])
      
      for(age in 2:(maxAge-1)){
        w[ind, age] <- (g[ind, age-1]/(lambdaInd[ind] -z[ind, age]))*w[ind, age-1] 
      }
      
      w[ind, maxAge] <- (g[ind, maxAge-1]/((lambdaInd[ind] -z[ind, maxAge])))*w[ind, maxAge-1] 
      
      #stable state distribution
      C[ind, 1:maxAge] <- w[ind, 1:maxAge]/sum(w[ind, 1:maxAge]) 
      #expected counts in each stage class
      for(age in 2:(maxAge)){
        ey[ind, age] <- (g[ind, age-1]*C[ind, age-1] + z[ind, age]*C[ind, age]) *selectivity[ind, age]
      }
      #ey[ind, maxAge] <- g[ind, maxAge-1]*C[ind, maxAge-1] + survivalProb[ind, maxAge]*C[ind, maxAge]
      ey[ind, 1] <- (lambdaInd[ind] - sum(ey[ind, 2:maxAge])) *selectivity[ind, 1]
      #ey[ind, 1] <- ( g[ind, 1]*C[ind, 1]) * fecundicity[ind, 1]
      
      # likelihood
       #theta[ind, 1:maxAge] <- rdirch(1,ey[ind, 1:maxAge])
      
       y[ind, 1:maxAge] <- rmultinom(1, prob = ey[ind, 1:maxAge], size = data$harvestCount[ind])
      
      #   # simulate new data
      #y_new[ind, 1:maxAge] ~ dmulti(ey[ind, 1:maxAge], harvestCount[ind])
      #   
      #for(age in 1:maxAge){
      #  g_obs[ind, age] <- (y[ind, age] + 0.5) * log(y[ind, age] + 0.5 / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
      #  g_new[ind, age] <- (y_new[ind, age] + 0.5) * log(y_new[ind, age] + 0.5 / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
      #  }
    }
  
  ret <- list(y = y, 
              beta = beta,
              lambdaInd = lambdaInd,
              #spawnProb = spawnProb,
              bsurvProbAge = bsurvProbAge)
  return(ret)
  
  }
  

simulatedDatasets <- lapply(seq(1,200,1), function(x){
  simData(data, constants, seed = x)
  })

save(simulatedDatasets, file = "result/simulatedDatasets.RData")

