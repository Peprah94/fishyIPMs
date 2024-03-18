# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs")


# Load formatted data
load("dataset/formattedDataList.RData")

#load results
load("modelFitting/fittedModel.RData")

# Load required packages
library(nimble)
library(dplyr)
library(readr)
library(reshape2)
library(ggplot2)

groupIndices <- dataList$ageAtHarvestData[, 1:3]

# Extract population growth rate
interestedVars <- rownames(fishModelMCMCrun$summary$all.chains)[grepl("lambda", rownames(fishModelMCMCrun$summary$all.chains))]

growthRate <- fishModelMCMCrun$summary$all.chains[interestedVars, ]%>%
  cbind(., groupIndices)

all <- growthRate%>%
  tidyr::expand(yearGrowthOcc, sex, lakeName)

growthRate <- growthRate %>% right_join(all)

indData <- dataList$indData

allGroups <- tidyr::expand(yearGrowthOcc = as.factor(unique(indData$yearGrowthOcc[!is.na(indData$yearGrowthOcc)])), 
                             sex = as.factor(unique(indData$sex[!is.na(indData$sex)])), 
                             lakeName = as.factor(unique(indData$lakeName[!is.na(indData$lakeName)])))
# Merge indData and all groups
growthRate <- growthRate%>%
  full_join(allGroups, 
             by = c("yearGrowthOcc", "sex", "lakeName"))

growthRate[is.na(growthRate)] <- 0


# Let's see the variation of the population growth rate
growthRate%>%
  ggplot(., aes(x = as.factor(yearGrowthOcc), y = Mean, group = as.factor(sex)))+
  geom_point(aes(x = as.factor(yearGrowthOcc), y = Mean, color = as.factor(sex)))+
  geom_line()+
  geom_ribbon(aes(ymin = `95%CI_low`, ymax = `95%CI_upp`, fill = as.factor(sex)), alpha = 0.2)+
  facet_wrap( ~ as.factor(lakeName))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  ylab("Population Growth Rate")+
  xlab("Year Growth Occured")
