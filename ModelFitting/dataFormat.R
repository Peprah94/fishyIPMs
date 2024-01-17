# Modelling the data

#load packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)

# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs")

# load Population data
popnData <- read_delim("dataset/Example1.csv", 
                       delim = ";", 
                       escape_double = FALSE, 
                       trim_ws = TRUE)

colnames(popnData)[c(2:5,8,10,12,16:25,29:30)] <- c("climateZone",
                           "lakeNumber",
                           "lineNumber",
                           "date",
                           "lakeName",
                           "gillnetNumber",
                           "gillnetMeshSize",
                           "speciesName",
                           "fishLength",
                           "fishWeight",
                           "fishFat",
                           "gonadWeight",
                           "sex",
                           "maturationStage",
                           "maturation",
                           "revisedMaturation",
                           "stomachContent",
                           "otoliothRadius",
                           "captureAge")

# remove NAs from some selected columns
popnData <- popnData[complete.cases(popnData[, c("maturationStage", "maturation")]),] %>%
  dplyr::mutate(maturationStage = as.numeric(maturationStage),
                maturation = as.numeric(maturation))

# Explore the relationship
r <- c(popnData[, 22])$maturationStage
s <- c(popnData[, 23])$maturation
table(s,r) #what is the relationship between maturation and maturation stage?


# Individual dataset
indData <- read_delim("dataset/Example2.csv", 
                       delim = ";", 
                       escape_double = FALSE, 
                       trim_ws = TRUE)

colnames(indData)[c(1:6,8:34)] <- c("lengthYrBefore",
                             "legthAtAgeThisYear",
                             "yearGrowthOcc",
                             "ageGrowthOcc",
                             "ageAtYr",
                             "fishID",
                             "climateZone",
                             "locationNumber",
                             "metersAboveSea",
                             "lineNumber",
                             "date",
                             "year",
                             "period",
                             "lakeName",
                             "fishAge",
                             "fishLength",
                             "fishWeight",
                             "sex",
                             "maturationStage",
                             "maturation",
                             "ageAtGrowth",
                             "growth",
                             "instanteneousGrowth",
                             "capturePerUnitEffort",
                             "temperatureMeanInMay",
                             "temperatureMeanInJune",
                             "temperatureMeanInJuly",
                             "precipitationSumInMay",
                             "precipitationSumInJune",
                             "precipitationSumInJuly",
                             "temperatureMeanInSummer",
                             "precipitationSumDuringSummer",
                             "NAOwinterIndex")

indData <- indData[complete.cases(indData[, c("maturationStage", "maturation")]),] %>%
  dplyr::mutate(maturationStage = as.numeric(maturationStage),
                maturation = as.numeric(maturation),
                presence = 1)


# catchment variables
catchmentVars <- read_csv("dataset/NO_all_vars_over_catchments.csv")
colnames(catchmentVars)  

lakeNames <- read_delim("dataset/lakes_used_for_model_prep.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)


# get age-at-harvest data
# goal is to get the number of individuals harvested at for each age, per year, sex and lake number
# consistent with Skelly et al (2023)

# data is collected in 2008
ageAtHarvestData <- dcast(indData, 
                          yearGrowthOcc + sex + lakeName ~ ageAtYr, 
                          value.var = "presence",
                          fun.aggregate = sum)%>%
  left_join(., indData, 
            by = c("yearGrowthOcc", "sex", "lakeName"),
            keep = FALSE,
            multiple = "first")



# Now I get the data to model sprawning
# Use the idea in Magnus et.al 2021
#Select females and use the maturation variable (0/1) to model the sprawning probability
index <- function(x) {ifelse(x > 0, 1, 0)}

sprawningAtHarvestData <- indData %>%
  dplyr::mutate(maturation = ifelse(sex == 1, 0, maturation))%>% # I assume 1 is males
  dcast(., 
        yearGrowthOcc + sex + lakeName ~ ageAtYr, 
        value.var = "maturation",
        fun.aggregate = sum)%>%
  left_join(., indData, 
            by = c("yearGrowthOcc", "sex", "lakeName"),
            keep = FALSE,
            multiple = "first")

# I intend to use the stock recruitment function.
# however, I don't have access to the number of eggs


# Catchment variables

#save data to use in fitting the models

dataList <- list(ageAtHarvestData = ageAtHarvestData,
                 sprawningData = sprawningAtHarvestData,
                 indData = indData,
                 popnData = popnData)
  
save(dataList, 
     file = "dataset/formattedDataList.RData")
