# Modelling the data

#load packages
library(nimble)
library(dplyr)
library(readr)
library(readxl)
library(reshape2)
library(tidyr)
library(RFishBC)
library(FSA)
library(tidyverse)

# set working directory.
#change this for your particular application
setwd("C:/GitHub/fishyIPMs/NewDataset")


# Import Individual data
indData <- read.delim("dataset/innlandsfisk.txt",
                      sep = ";")

# select the species ID of interest from the data
focal_speciesid <- 5
indData <- indData[indData$FK_ArtID==focal_speciesid,]


#remove NAs for Length and Radius 
indData <- indData[complete.cases(indData$Radius), ]
indData <- indData[complete.cases(indData$Lengde), ]

#Check fish data for different lakes
table(indData$FK_InnsjoNr)

# Some lakes have very few individuals across time(less than 10), 
#remove these lakes as they dont seem representative for the species
# Remove rows with less than min_rows_per_category per category
indData <- indData %>%
  group_by(FK_InnsjoNr) %>%
  filter(n() >= 9) %>%
  ungroup()

plot(indData$Lengde ~ indData$Radius)

#Seem to be two types of measurements, just remove the "outliers" for now
indData <- indData[indData$Radius<400,]

#Cheking for strange length Radius correlations among lakes and years
data_check <- indData %>%
  group_by(FK_InnsjoNr, Aar) %>%
  summarize(R_squared = summary(lm(Lengde ~ Radius))$r.squared,
            Intercept = coef(lm(Lengde ~ Radius))[1])
print(data_check)

#decide on cut-off on 0.8
ok_lakes_years <- data_check[data_check$R_squared>0.8, ]

#Remove r_squared of 1, which seems wronge
ok_lakes_years <- ok_lakes_years[!ok_lakes_years$R_squared==1, ]

indData <- indData[indData$FK_InnsjoNr %in% ok_lakes_years$FK_InnsjoNr & indData$Aar %in% ok_lakes_years$Aar, ]
#plot(backcalc_data$Lengde~backcalc_data$Radius)

#setting up data for RFishBC (renaming variables)
indData <- indData %>%
  rename(id = InnlandsfiskID,
         agecap = AlderSkjell,
         radcap = Radius,
         rad1 = Radius1,
         rad2 = Radius2,
         rad3 = Radius3,
         rad4 = Radius4,
         rad5 = Radius5,
         rad6 = Radius6,
         rad7 = Radius7,
         rad8 = Radius8,
         rad9 = Radius9,
         rad10 = Radius10,
         rad11 = Radius11,
         rad12 = Radius12,
         rad13 = Radius13,
         rad14 = Radius14,
         rad15 = Radius15,
         rad16 = Radius16,
         rad17 = Radius17,
         rad18 = Radius18,
         rad19 = Radius19,
         rad20 = Radius20)

##Back calculating data using fraser lee per lake 
indData <- indData %>%
  group_by(FK_InnsjoNr) %>%
  do(backCalc(., Lengde, BCM = "FRALE",
              inFormat = "wide", outFormat = "long", digits = 0)) %>%
  ungroup()

# Print the resulting dataframe
print(indData)

###
#Rename variables
indData <- indData %>%
  rename(LengdeValder = bclen,
         alder_aar = ann)

#including year of length increment (specific year affecting specific growth)
indData$vekstaar <- indData$Aar-((indData$agecap - indData$alder_aar+1))


#select variables of interest and change column names
indData <- indData[, c(1, 4, 6:11, 21, 40:42)]
colnames(indData) <- c("fishID", 
                        "InnsjoNr",
                       "date",
                       "year",
                       "month",
                       "day",
                       "runNr",
                       "speciesID",
                       "sex",
                       "ageAtYear",
                       "lengthAtYear",
                       "yearGrowthOccured")

# import data with catchment variables
catchmentVars <- read_csv("dataset/timeseries_combined_vars.csv")%>%
                  dplyr::mutate(innsjo_nr = as.character(innsjo_nr))%>%
                  dplyr::select(-c("vassdragNr"))%>% 
                  distinct(innsjo_nr, year, .keep_all = TRUE)%>%
                   rowwise() %>%
                  mutate(summerTemperature = mean(c_across(c('may_temp_mean','jun_temp_mean', 'jul_temp_mean', 'aug_temp_mean')), na.rm = TRUE),
                         summerPrecipitation = mean(c_across(c('may_total_precip','jun_total_precip', 'jul_total_precip', 'aug_total_precip')), na.rm = TRUE),
                         summerSnowDepth = mean(c_across(c('may_snow_depth','jun_snow_depth', 'jul_snow_depth', 'aug_snow_depth')), na.rm = TRUE))
colnames(catchmentVars)[1] <- c("InnsjoNr") # ensure that the lakenr are the same for all datasets

# import lake IDs
lakeInfo <- read_excel("dataset/Oversikt-lokaliteter-lagt inn i basen.xls")%>%
            dplyr::select("Loknr", "InnsjoNr", "Loknavn") 


# Add lake Nr and name to the individual fish dataset
# Catchment dataset has been formatted to the regions and years we need.
indDataWithCatchments <- indData %>%
  filter(sex <= 2)%>%
  filter(InnsjoNr %in% catchmentVars$InnsjoNr)%>%
  filter(year %in% catchmentVars$year)%>%
  inner_join(catchmentVars, by = c("InnsjoNr", "year"), keep = FALSE)


# Put everything together
#Estimate age at year, year growth occured
# indDataWithCatchments <- indData1%>%
#   select(c(InnsjoNr, year, fishID, starts_with("radius")))%>% 
#   reshape2::melt(., id.vars = c("InnsjoNr", "year", "fishID"), 
#                  variable.name = "radius",
#                  value.name = "lengthAtAge")%>%
#   dplyr::group_by(fishID)%>% # we group the data by the fishID
#   dplyr::mutate(maxAge = sum(!is.na(lengthAtAge)))%>% # maxAge is the sum of radius values that are not NAs
#   #note that we assume radius represents the length at age
#   dplyr::filter(., !is.na(lengthAtAge))%>% # we take out the NAs before calculating the age at year
#   dplyr::mutate(yearGrowthOccured = ifelse(maxAge > 1, year + 1- seq(0:(maxAge-1)), year), #year growth occured is backdating the year of harvest with the number of columns with radius observations
#                 ageAtYear = ifelse(maxAge > 1, maxAge  - seq(0:(maxAge-1)), maxAge-1)) %>%
#   ungroup()%>%
#   left_join(., indData1, by = c("InnsjoNr", "year", "fishID"))




# get age-at-harvest data
# goal is to get the number of individuals harvested at for each age, per year, sex and lake number
# consistent with Skelly et al (2023)

# data is collected in 2008
# create the age at harvest data
indDataWithCatchments$presence <- 1
ageAtHarvestData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("ageAtYear", "sex", "InnsjoNr", "speciesID")]),]%>%
  dcast(., 
        year + sex + InnsjoNr  ~ ageAtYear, 
        value.var = "presence",
        fun.aggregate = sum)%>%
  left_join(., indDataWithCatchments, 
            by = c("year", "sex", "InnsjoNr"),
            keep = FALSE,
            multiple = "first")



# Now I get the data to model sprawning
# Use the idea in Magnus et.al 2021
#Select females and use the maturation variable (0/1) to model the sprawning probability
index <- function(x) {ifelse(x > 0, 1, 0)}

sprawningAtHarvestData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("ageAtYear", "sex", "InnsjoNr", "speciesID")]),]%>%
  dplyr::mutate(maturation = ifelse(ageAtYear > 2, 0,1))%>% # I assume maturation occurs when age is grater than 1
  dcast(., 
        year + sex + InnsjoNr ~ ageAtYear, 
        value.var = "maturation",
        fun.aggregate = sum)%>%
  left_join(., indDataWithCatchments, 
            by = c("year", "sex", "InnsjoNr"),
            keep = FALSE,
            multiple = "first")

# I intend to use the stock recruitment function.
# however, I don't have access to the number of eggs


# Length at age dataset
lengthAtAgeData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("speciesID", "sex", "InnsjoNr", "ageAtYear")]),]%>%
  #dplyr::mutate(maturation = ifelse(sex == 1, 0, maturation))%>% # I assume 1 is males
  dcast(., 
        year + sex + InnsjoNr ~ ageAtYear, 
        value.var = "lengthAtYear",
        fun.aggregate = mean)%>%
  left_join(., indDataWithCatchments, 
            by = c("year", "sex", "InnsjoNr"),
            keep = FALSE,
            multiple = "first")
#change the NANs into NAs
lengthAtAgeData[is.na(lengthAtAgeData)] <- NA

#save data to use in fitting the models

dataList <- list(ageAtHarvestData = ageAtHarvestData,
                 sprawningData = sprawningAtHarvestData,
                 indData = indDataWithCatchments,
                 #popnData = popnData,
                 lengthAtAgeData = lengthAtAgeData)

save(dataList, 
     file = "dataset/formattedDataList.RData")

# check the catchment variables and correlation with other variables

#lengthAtAgeThisYear
#GGally::ggpairs(ageAtHarvestData[, c(1:3, 15, 36:45)])

# Catchment variables
#GGally::ggpairs(ageAtHarvestData[, c(53:68)])


