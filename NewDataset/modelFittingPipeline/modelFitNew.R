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
  
  for(i in 1:13){
    bsurv[i] ~dnorm(0, sd = 10)
  }
  
  sdbsurvAge ~ dunif(0.01, 10)
  sdbsurvLake ~ dunif(0.01, 10)
  sdbsurvSex ~ dunif(0.01, 10)
  
  for(i in 1:maxAge){
    bsurvAge[i]~dnorm(0, sd = sdbsurvAge)
  }
  
  # variable selection parameters
  for(i in 1:12){
    psi[i] ~ dunif(0,1)
  }
  
  #lake effect
  for(i in 1:nLakes){
    bsurvLake[i] ~dnorm(0, sd = sdbsurvLake)
  }
  
  
  for(i in 1:nSex){
    bsurvSex[i] ~ dnorm(0, sd = sdbsurvSex)
  }
  
  for(ind in 1:nInds){
    for(age in 1:maxAge){ #loop through the age
      #logit(survivalProb[ind, age]) <- bsurv[1] + bsurvAge[age] + bsurvLake[lakeID[ind]] + bsurvSex[sexID[ind]] + psi[1] * bsurv[2]* year[ind]  + psi[2] * bsurv[3]*forest[ind] + psi[3] * bsurv[4]*pastures[ind] + psi[4] * bsurv[5]*popnDensity[ind] + psi[5] * bsurv[6]*summerTemp[ind] + psi[6] * bsurv[7]*summerPrec[ind] + psi[7] * bsurv[8]*summerSnowDepth[ind]  + psi[8] * bsurv[9]*forestsq[ind] + psi[9] * bsurv[10]*pasturessq[ind] + psi[10] * bsurv[11]*popnDensitysq[ind] + psi[11] * bsurv[12]*summerTempsq[ind] + psi[12] * bsurv[13]*summerPrecsq[ind] + psi[13] * bsurv[14]*summerSnowDepthsq[ind] 
      logit(survivalProb[ind, age]) <- bsurv[1] + bsurvAge[age] + bsurvLake[lakeID[ind]] + bsurvSex[sexID[ind]] + psi[1] * bsurv[2]* year[ind]  + psi[2] * bsurv[3]*forest[ind] + psi[3] *bsurv[4]*pastures[ind] + psi[4] *bsurv[5]*popnDensity[ind] + psi[5] *bsurv[6]*summerTemp[ind] + psi[6] *bsurv[7]*summerSnowDepth[ind]  + psi[7] *bsurv[8]*forestsq[ind] + psi[8] *bsurv[9]*pasturessq[ind] + psi[9] *bsurv[10]*popnDensitysq[ind] + psi[10] *bsurv[11]*summerTempsq[ind] + psi[11] *bsurv[12]*summerSnowDepthsq[ind] + psi[12] *bsurv[13]*WPUE[ind]
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
   # Amat[maxAge, maxAge, ind] <- survivalProb[ind, maxAge]

  }
  #} 
  
  # Growth rate (lambda) for each individual
  for(ind in 1:nInds){
    #lambda[ind] <- exp(mean(r.annual[u:T, ind])) 
    
    # Asymptotic population growth rate
    lambda[ind] <- nimbleLambdaEstimation(Amat[1:maxAge, 1:maxAge, ind])
    # Proportion of stable distribution
    # check skelly et. al (2023) page 1066 for equations
    
    # w = proportion to stable age
    # z = probability of remaining in the same age class
    # g = probability of surviving to the next stage class
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
    w[ind, maxAge] <- (g[ind, maxAge-1]/(lambda[ind] -z[ind, maxAge]))*w[ind, maxAge-1] 
    
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
    
    #   # simulate new data
      #y_new[ind, 1:maxAge] ~ dmulti(ey[ind, 1:maxAge], harvestCount[ind])
    #   
      #for(age in 1:maxAge){
      #  g_obs[ind, age] <- (y[ind, age] + 0.5) * log(y[ind, age] + 0.5 / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
      #  g_new[ind, age] <- (y_new[ind, age] + 0.5) * log(y_new[ind, age] + 0.5 / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
      #  }
  }
  
  # Model checking
  # for(ind in 1:nInds){
  #   # simulate new data
  #   y_new[ind, 1:maxAge] ~ dmulti(theta[ind, 1:maxAge], harvestCount[ind])
  #   
  #   for(age in 1:maxAge){
  #     g_obs[ind, age] <- y[ind, 1:age] * log(y[ind, 1:age] / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
  #     g_new[ind, age] <- y_new[ind, 1:age] * log(y_new[ind, 1:age] / (harvestCount[ind] * ey[ind, age] / sum(ey[ind, 1:maxAge])))
  #     }
  # }
  # 
   # compute statistc
    #cs_obs <- 2 * sum(g_obs[1:nInds, 1:maxAge])
    #cs_new <- 2 * sum(g_new[1:nInds, 1:maxAge])
    #bayesPval <- step(cs_new - cs_obs )
  # 
  # Average survival probability
  
})

stagedData <- rowSums(as.matrix(dataList$ageAtHarvestData[ ,6:17])) %>%
  cbind(as.matrix(dataList$ageAtHarvestData[ ,4:5]), .)

wpue = as.numeric(scale(dataList$ageAtHarvestData$WPUE))
wpue[is.na(wpue)] <- mean(wpue, na.rm = TRUE)
# Data for model

data = list(
 # ageAtYear = as.numeric(dataList$sprawningData$ageAtYear),
  lengthAtAgeThisYear = lengthAtAgePosteriorSummary$data,
  year = as.numeric(scale(dataList$ageAtHarvestData$year)),
  y = as.matrix(stagedData),
  harvestCount = rowSums(as.matrix(dataList$ageAtHarvestData[4:17])),
 forest = as.numeric(scale(dataList$ageAtHarvestData$`Coniferous forest`)),
pastures = as.numeric(scale(dataList$ageAtHarvestData$Pastures)),
 summerTemp = as.numeric(scale(dataList$ageAtHarvestData$meanSummerTemp)),
 #summerPrec = as.numeric(scale(dataList$ageAtHarvestData$summerPrecipitation)),
 summerSnowDepth = as.numeric(scale(dataList$ageAtHarvestData$meanWinterSnow)),
 popnDensity = as.numeric(scale(dataList$ageAtHarvestData$pop_density)),
 forestsq = as.numeric((scale(dataList$ageAtHarvestData$`Coniferous forest`)))^2,
 pasturessq = as.numeric((scale(dataList$ageAtHarvestData$Pastures)))^2,
 summerTempsq = as.numeric((scale(dataList$ageAtHarvestData$meanSummerTemp)))^2,
 #summerPrecsq = as.numeric((scale(dataList$ageAtHarvestData$summerPrecipitation)))^2,
 summerSnowDepthsq = as.numeric((scale(dataList$ageAtHarvestData$meanWinterSnow)))^2,
 popnDensitysq = as.numeric((scale(dataList$ageAtHarvestData$pop_density)))^2,
 WPUE = wpue
)

constants = list(
  nInds = nrow(dataList$ageAtHarvestData),
  sexID = as.numeric(dataList$sprawningData$sex),
  nSex = length(unique(as.numeric(dataList$sprawningData$sex))),
  lakeID = as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)),
  nLakes = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$InnsjoNr)))),
  maxAge = dim(data$y)[2]#,
  #nSpecies = length(unique(as.numeric(as.factor(dataList$ageAtHarvestData$speciesID))))
)

inits <- list(
  psi = rep(0.5, 12),
  bsurv = rnorm(13, 0, 1),
  bsurvSex = rnorm(constants$nSex, 0, 1),
  bsurvAge = rnorm(constants$maxAge, 0, 1),
  bsurvLake = rnorm(constants$nLakes, 0,1),
  sdbsurvSex = 1,
  sdbsurvAge = 1,
  sdbsurvLake = 1,
  Amat = array(0, dim = c(constants$maxAge,constants$maxAge, constants$nInds)),
  spawnProb = 0.1,
  theta = matrix((rep(1, constants$maxAge)/constants$maxAge), nrow = constants$nInds, ncol = constants$maxAge)
)



# nimbleModel
fishModel <- nimbleModel(code,
                         data = data,
                         constants = constants,
                         inits = inits)

# compile nimbleModel
fishModelCompiled <- compileNimble(fishModel)

# check if any lambda is 0 or NA
which(fishModelCompiled$lambda == 0)
which(is.na(fishModelCompiled$lambda))

# Configure the model
fishModelConfigured <- configureMCMC(fishModelCompiled,
                                     monitors = c("bsurv", 
                                                  "bsurvLake", 
                                                  "bsurvAge", 
                                                  "bsurvSex",
                                                  "sdbsurvLake", 
                                                  "sdbsurvAge",  
                                                  "sdbsurvSex",
                                                  "lambda", 
                                                  "survivalProb",
                                                  "spawnProb", 
                                                  "psi"#,
                                                  #"bayesPval"
                                                  )
)

# Build and run model
fishModelBuilt <- buildMCMC(fishModelConfigured)

# compile built model
fishModelCompiled2 <- compileNimble(fishModelBuilt,
                                    project = fishModelCompiled)


#run MCMC
fishModelMCMCrun <- runMCMC(fishModelCompiled2,
                            niter = 10000,
                            nchains = 2,
                            nburnin = 5000,
                            thin = 2,
                            setSeed = TRUE,
                            samples = TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE) 

#check summary
fishModelMCMCrun$summary$all.chains

save(fishModelMCMCrun, file = "result/fittedModelWithOverdispersion.RData")

# currently the model runs. Will have to think about the model structure in details.