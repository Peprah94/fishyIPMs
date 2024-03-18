# Modelling the data

#load packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)
library(tidyr)
# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs")
# catchment variables
catchmentVars <- read_csv("dataset/model_catchment_vars.csv")
colnames(catchmentVars)[2] <- "lakeID"  



# Individual dataset
indData <- read_csv("dataset/Example2_vatnLnr.csv")

colnames(indData)[c(1:6,8:35)] <- c("lengthYrBefore",
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
                                    "lakeID",
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


#merge indData and catchment variables

indDataWithCatchments <- indData %>%
  right_join(catchmentVars, by = c("lakeID", "year"))%>%
  dplyr::mutate(maturationStage = as.numeric(maturationStage),
                maturation = as.numeric(maturation),
                presence = 1)

#unique(indDataWithCatchments$year)
#unique(indData$year)

#unique(indDataWithCatchments$yearGrowthOcc[!is.na(indDataWithCatchments$yearGrowthOcc)])
#unique(indDataWithCatchments$sex)
#unique(indDataWithCatchments$lakeName)







 table(indDataWithCatchments$yearGrowthOcc, indDataWithCatchments$lakeName)
  
#dcast(indDataWithCatchments, formula = yearGrowthOcc + sex + lakeName ~ ageAtYr, 
#      value.var = "presence", 
 #     drop = FALSE,
 #     fill = 0)





# lakeNames <- read_delim("dataset/lakes_used_for_model_prep.csv", 
#                         delim = ";", 
#                         escape_double = FALSE, 
#                         trim_ws = TRUE)


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






# get age-at-harvest data
# goal is to get the number of individuals harvested at for each age, per year, sex and lake number
# consistent with Skelly et al (2023)

# data is collected in 2008
# create the age at harvest data
ageAtHarvestData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("ageAtYr", "sex", "lakeName")]),]%>%
                  
  
  dcast(., 
                          yearGrowthOcc + sex + lakeName ~ ageAtYr, 
                          value.var = "presence",
                          fun.aggregate = sum)%>%
  left_join(., indDataWithCatchments, 
            by = c("yearGrowthOcc", "sex", "lakeName"),
            keep = FALSE,
            multiple = "first")



# Now I get the data to model sprawning
# Use the idea in Magnus et.al 2021
#Select females and use the maturation variable (0/1) to model the sprawning probability
index <- function(x) {ifelse(x > 0, 1, 0)}

sprawningAtHarvestData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("ageAtYr", "sex")]),]%>%
  #dplyr::mutate(maturation = ifelse(sex == 1, 0, maturation))%>% # I assume 1 is males
  dcast(., 
        yearGrowthOcc + sex + lakeName ~ ageAtYr, 
        value.var = "maturation",
        fun.aggregate = sum)%>%
  left_join(., indDataWithCatchments, 
            by = c("yearGrowthOcc", "sex", "lakeName"),
            keep = FALSE,
            multiple = "first")

# I intend to use the stock recruitment function.
# however, I don't have access to the number of eggs


# Catchment variables

#save data to use in fitting the models

dataList <- list(ageAtHarvestData = ageAtHarvestData,
                 sprawningData = sprawningAtHarvestData,
                 indData = indDataWithCatchments,
                 popnData = popnData)
  
save(dataList, 
     file = "dataset/formattedDataList.RData")

# check the catchment variables and correlation with other variables

#lengthAtAgeThisYear
GGally::ggpairs(ageAtHarvestData[, c(1:3, 15, 36:45)])

# Catchment variables
GGally::ggpairs(ageAtHarvestData[, c(53:68)])


