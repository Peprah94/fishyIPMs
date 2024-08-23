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
#setwd("C:/GitHub/fishyIPMs/NewDataset")


# Import Individual data
indData <- read.delim("dataset/innlandsfisk.txt",
                      sep = ";")

# select the species ID of interest from the data
focal_speciesid <- 5
indData <- indData[indData$FK_ArtID==focal_speciesid,]

#Incklude biotic variables
#Include information whether a fish population is sympatric or allopatric, 
# and also indicate if large predators are present (i.e., pike)

allopatric_species <- indData %>%
  group_by(FK_InnsjoNr) %>%
  summarise(num_species = n_distinct(FK_ArtID)) %>%
  filter(num_species == 1) %>%
  select(FK_InnsjoNr)

indData$is_allopatric <- ifelse(indData$FK_InnsjoNr %in% allopatric_species$FK_InnsjoNr, "allopatric", "sympatric")

# add pike info
indData <- indData %>%
  group_by(FK_InnsjoNr) %>%
  mutate(pike_presence = if_else("37" %in% unique(FK_ArtID), "Yes", "No"))

#Calculate wpue per fishing day, assuming same effort across lakes and year (one night per event)

indData <- indData %>%
  group_by(FK_InnsjoNr, Aar, Maaned,Dag) %>%
  mutate(WPUE_day = sum(Vekt)) %>%
  ungroup()

#average WPUE over a year if multiple fishing days
indData <- indData %>%
  group_by(FK_InnsjoNr, Aar) %>%
  mutate(WPUE_avg=mean(WPUE_day))%>%
  ungroup()

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
indData <- indData[, c(1, 4, 6:10, 21, 40:41, 43, 44:46)]
colnames(indData) <- c("fishID", 
                        "InnsjoNr",
                       "date",
                       "year",
                       "month",
                       "day",
                       "runNr",
                       "sex",
                       "allopatric",
                       "pikePresence",
                       "WPUE",
                       "ageAtYear",
                       "lengthAtYear",
                       "yearGrowthOccured")


#Select the data without missing ageAtYear, sex and lakeID, and sex either 1 or 2
indData <- indData[complete.cases(indData[, c("ageAtYear", "sex", "InnsjoNr")]),]%>%
  #filter(year < 2010) %>%
  filter(sex < 3)

length(unique(indData$fishID))
# Two years do not have any WPUE value, so I will take the mean of the 1991 and 1993 values for those
indData$WPUE[which(is.na(indData$WPUE))] <- mean(indData$WPUE[indData$year %in% c("1991", "1993")])


# Create covariates for future predictions
yearPred <- seq(2013, 2030, 1) 
sexPred <- c(1,2)
InnsjoNrPred <- unique(indData$InnsjoNr)[!is.na(unique(indData$InnsjoNr))]
ageAtYearPred <- seq(1,14, 1)

indDataPred <- tidyr::crossing(yearPred, 
                               sexPred, 
                               InnsjoNrPred, 
                               ageAtYearPred)

colnames(indDataPred) <- c("year", 
                           "sex",
                           "InnsjoNr",
                           "ageAtYear")

indDataPred$month <- NA
indDataPred$day <- NA
indDataPred$runNr <- NA
indDataPred$allopatric <- "allopatric"
indDataPred$pikePresence <- "No"
indDataPred$fishID <- NA
indDataPred$yearGrowthOccured <- NA
indDataPred$lengthAtYear <- NA
indDataPred$WPUE <- NA
indDataPred$date <- NA

indData <- rbind(indData,
                 indDataPred)


# import data with climatic variables

climaticVars <- read_csv("dataset/expanded_climate_vars.csv")
names(climaticVars)

climaticVars$meanTempSummer <-rowMeans(climaticVars[c("may_temp_mean",
                                                      "jun_temp_mean",
                                                      "jul_temp_mean",
                                                      "aug_temp_mean")])

climaticVars$snowDecToMarch <-rowMeans(climaticVars[c("dec_snow_depth", 
                                                      "jan_snow_depth", 
                                                      "feb_snow_depth", 
                                                      "mar_snow_depth")])

climaticVars <- climaticVars%>%
  select(c("vassdragNr", "innsjo_nr", "year", "meanTempSummer", "snowDecToMarch"))%>%
  rename(InnsjoNr = innsjo_nr)



##match based on lakeID and year
#mean summer temp
indData$meanSummerTemp <- climaticVars$meanTempSummer[match(paste(indData$InnsjoNr, indData$year), 
                                                            paste(climaticVars$InnsjoNr,climaticVars$year))]

#Snow in spring and autumn 
indData$meanWinterSnow <- climaticVars$snowDecToMarch[match(paste(indData$InnsjoNr, indData$year), 
                                                            paste(climaticVars$InnsjoNr, climaticVars$year))]

#Format Catchment Variables
catchmentVars <- read_csv("dataset/timeseries_combined_vars.csv")%>%
                  dplyr::mutate(innsjo_nr = as.character(innsjo_nr))

indData$forest <- catchmentVars$`Coniferous forest`[match(paste(indData$InnsjoNr, indData$year),
                                                            paste(catchmentVars$innsjo_nr, catchmentVars$year))]

indData$popnDensity <- catchmentVars$pop_density[match(paste(indData$InnsjoNr, indData$year),
                                                               paste(catchmentVars$innsjo_nr, catchmentVars$year))]

#%>%
                 # dplyr::select(-c("vassdragNr"))%>% 
                 # distinct(innsjo_nr, year, .keep_all = TRUE)%>%
                 # group_by(year)%>%
                 # summarise_all(., .funs = "mean")%>%
                 # ungroup()%>%
                 # select(-c("innsjo_nr"))

#colnames(catchmentVars)[1] <- c("InnsjoNr") # ensure that the lakenr are the same for all datasets

# import lake IDs
lakeInfo <- read_excel("dataset/Oversikt-lokaliteter-lagt inn i basen.xls")%>%
            dplyr::select("Loknr", "InnsjoNr", "Loknavn", "Y") 


# Add lake Nr and name to the individual fish dataset
# Catchment dataset has been formatted to the regions and years we need.
# indDataWithCatchments <- indData %>%
#   filter(sex <= 2)%>%
#  # filter(InnsjoNr %in% catchmentVars$InnsjoNr)%>%
#   filter(year %in% catchmentVars$year)%>%
#   inner_join(catchmentVars, by = c("year"), keep = FALSE)
#"InnsjoNr",

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
indDataWithCatchments <- indData
indDataWithCatchments$presence <- 1
ageAtHarvestData <- indDataWithCatchments %>%
  dcast(., 
        year + sex + InnsjoNr  ~ ageAtYear, 
        value.var = "presence",
        fun.aggregate = sum)%>%
  left_join(., indDataWithCatchments, 
            by = c("year", "sex", "InnsjoNr"),
            keep = FALSE,
            relationship = "many-to-many",
            multiple = "all")%>%
  group_by(year, sex, InnsjoNr)%>%
  summarise_all(., .funs = "mean") %>%
  left_join(., lakeInfo, by = "InnsjoNr")

# Set the latitude of the covariate
ageAtHarvestData[ageAtHarvestData$InnsjoNr %in% "1649",]$Y <- 6810000
ageAtHarvestData[ageAtHarvestData$InnsjoNr %in% "29590",]$Y <- 6842209

# Remove the data without latitude
ageAtHarvestData <- ageAtHarvestData[complete.cases(ageAtHarvestData[, c("Y")]), ]


# Simulate covariates for the predictions
ageAtHarvestDataObs <- ageAtHarvestData %>%
  filter(year < 2013) # We have data up to 2012

ageAtHarvestDataObs <- ageAtHarvestDataObs[complete.cases(ageAtHarvestDataObs[, c("meanSummerTemp", "meanWinterSnow", "Y")]), ]

ageAtHarvestDataPred <- ageAtHarvestData %>%
  filter(year > 2012) # We want to predict from 2013 to 2030

InnsjoNrPred <- InnsjoNrPred[InnsjoNrPred %in% ageAtHarvestDataObs$InnsjoNr]

  for(l in seq_along(InnsjoNrPred)){
    print(l)
    indx <- which(ageAtHarvestDataPred$InnsjoNr %in% InnsjoNrPred[l])
    tp <- which(ageAtHarvestDataObs$InnsjoNr %in% InnsjoNrPred[l])
 print(length(tp))
    
    if(length(indx) > 1){
    ageAtHarvestDataPred[indx, "meanSummerTemp"] <- rnorm(length(indx),
                                            mean = mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp), na.rm = TRUE),
                                            sd = sd(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp),na.rm = TRUE))
  
    ageAtHarvestDataPred[indx, "meanWinterSnow"] <- rnorm(length(indx),
                                                          mean = mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE),
                                                          sd = sd(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE))
    
    ageAtHarvestDataPred[indx, "forest"] <- rnorm(length(indx),
                                                          mean = mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE),
                                                          sd = sd(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE))
    
    
    ageAtHarvestDataPred[indx, "popnDensity"] <- rnorm(length(indx),
                                                          mean = mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE),
                                                          sd = sd(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE))
   
     } else {
       ageAtHarvestDataPred[indx, "meanSummerTemp"] <- mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp), na.rm = TRUE)
       ageAtHarvestDataPred[indx, "meanWinterSnow"] <-  mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE)
       ageAtHarvestDataPred[indx, "forest"] <- mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE)
       ageAtHarvestDataPred[indx, "popnDensity"] <- mean(c(ageAtHarvestDataObs[ageAtHarvestDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE)
      
    }
    
   
     }



# Fill the NAs for WPUE, forest and population density
fillNAwithMean <- function(x, data){ # x is the colum index for the covariate of interest 
  na_index <-  which(is.na(data[, x] ))  
  yrs <- data[na_index, ]$year
  for(i in seq_along(yrs)){
    filterData <- data%>%
                  filter(year == yrs[i])
    data[na_index[i],x] <- mean(c(filterData[,x][[1]]), na.rm = TRUE)
  }
  return(data)
}

ageAtHarvestData <- rbind(ageAtHarvestDataObs,
                          ageAtHarvestDataPred)%>% 
  mutate_all(~ifelse(is.nan(.), NA, .))
# Fill NAs for forest
ageAtHarvestData <- fillNAwithMean(31, ageAtHarvestData)
# Fill NAs for population density
ageAtHarvestData <- fillNAwithMean(32, ageAtHarvestData)
# Fill NAs for WPUE
ageAtHarvestData <- fillNAwithMean(25, ageAtHarvestData)
ageAtHarvestData$WPUE[is.na(ageAtHarvestData$WPUE)] <- mean(ageAtHarvestData$WPUE[!is.na(ageAtHarvestData$WPUE)])
# Check if they are indeed filled
which(is.na(ageAtHarvestData$forest))
which(is.na(ageAtHarvestData$popnDensity))
which(is.na(ageAtHarvestData$WPUE))


# Put the ageAtHarvest data together

# ageAtHarvestData <- ageAtHarvestData%>%
#   mutate_at(., .vars = c("meanSummerTemp", "meanWinterSnow", "forest", "popnDensity"), .funs = "scale")%>%
#   ungroup()

# Now I get the data to model sprawning
# Use the idea in Magnus et.al 2021
#Select females and use the maturation variable (0/1) to model the sprawning probability
index <- function(x) {ifelse(x > 0, 1, 0)}

sprawningAtHarvestData <- indDataWithCatchments[complete.cases(indDataWithCatchments[, c("ageAtYear", "sex", "InnsjoNr")]),]%>%
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
lengthAtAgeData <- indDataWithCatchments %>%
  #dplyr::mutate(maturation = ifelse(sex == 1, 0, maturation))%>% # I assume 1 is males
  dcast(., 
        year + sex + InnsjoNr ~ ageAtYear, 
        value.var = "lengthAtYear",
        fun.aggregate = mean)%>%
  left_join(., indDataWithCatchments, 
            by = c("year", "sex", "InnsjoNr"),
            keep = FALSE,
            relationship = "many-to-many",
            multiple = "all")%>%
  group_by(year, sex, InnsjoNr)%>%
  summarise_all(., .funs = "mean") %>%
  left_join(., lakeInfo, by = "InnsjoNr")

#change the NANs into NAs
lengthAtAgeData[is.na(lengthAtAgeData)] <- NA

# Set the latitude of the covariate
lengthAtAgeData[lengthAtAgeData$InnsjoNr %in% "1649",]$Y <- 6810000
lengthAtAgeData[lengthAtAgeData$InnsjoNr %in% "29590",]$Y <- 6842209

# Remove the data without latitude
lengthAtAgeData <- lengthAtAgeData[complete.cases(lengthAtAgeData[, c("Y")]), ]


# Simulate covariates for the predictions
lengthAtAgeDataObs <- lengthAtAgeData %>%
  filter(year < 2013) # We have data up to 2012

lengthAtAgeDataObs <- lengthAtAgeDataObs[complete.cases(lengthAtAgeDataObs[, c("meanSummerTemp", "meanWinterSnow", "Y")]), ]

lengthAtAgeDataPred <- lengthAtAgeData %>%
  filter(year > 2012) # We want to predict from 2013 to 2030

InnsjoNrPred <- InnsjoNrPred[InnsjoNrPred %in% lengthAtAgeDataObs$InnsjoNr]

for(l in seq_along(InnsjoNrPred)){
  print(l)
  indx <- which(lengthAtAgeDataPred$InnsjoNr %in% InnsjoNrPred[l])
  tp <- which(lengthAtAgeDataObs$InnsjoNr %in% InnsjoNrPred[l])
  print(length(tp))
  
  if(length(indx) > 1){
    lengthAtAgeDataPred[indx, "meanSummerTemp"] <- rnorm(length(indx),
                                                          mean = mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp), na.rm = TRUE),
                                                          sd = sd(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp),na.rm = TRUE))
    
    lengthAtAgeDataPred[indx, "meanWinterSnow"] <- rnorm(length(indx),
                                                          mean = mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE),
                                                          sd = sd(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE))
    
    lengthAtAgeDataPred[indx, "forest"] <- rnorm(length(indx),
                                                  mean = mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE),
                                                  sd = sd(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE))
    
    
    lengthAtAgeDataPred[indx, "popnDensity"] <- rnorm(length(indx),
                                                       mean = mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE),
                                                       sd = sd(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE))
    
  } else {
    lengthAtAgeDataPred[indx, "meanSummerTemp"] <- mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanSummerTemp" ]$meanSummerTemp), na.rm = TRUE)
    lengthAtAgeDataPred[indx, "meanWinterSnow"] <-  mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "meanWinterSnow" ]$meanWinterSnow),na.rm = TRUE)
    lengthAtAgeDataPred[indx, "forest"] <- mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "forest" ]$forest),na.rm = TRUE)
    lengthAtAgeDataPred[indx, "popnDensity"] <- mean(c(lengthAtAgeDataObs[lengthAtAgeDataObs$InnsjoNr%in% InnsjoNrPred[l], "popnDensity" ]$popnDensity),na.rm = TRUE)
    
  }
  
  
}


lengthAtHarvestData <- rbind(lengthAtAgeDataObs,
                             lengthAtAgeDataPred)%>% 
  mutate_all(~ifelse(is.nan(.), NA, .))
# Fill NAs for forest
lengthAtHarvestData <- fillNAwithMean(31, lengthAtHarvestData)

# Fill NAs for population density
lengthAtHarvestData <- fillNAwithMean(32, lengthAtHarvestData)

# Fill NAs for WPUE
lengthAtHarvestData <- fillNAwithMean(25, lengthAtHarvestData)
lengthAtHarvestData$WPUE[is.na(lengthAtHarvestData$WPUE)] <- mean(lengthAtHarvestData$WPUE[!is.na(lengthAtHarvestData$WPUE)])
# Check if they are indeed filled
which(is.na(lengthAtAgeDataPred$forest))
which(is.na(lengthAtAgeDataPred$popnDensity))
which(is.na(lengthAtAgeDataPred$WPUE))

# lengthAtHarvestData <- lengthAtHarvestData%>%
#   mutate_at(., .vars = c("meanSummerTemp", "meanWinterSnow", "forest", "popnDensity"), .funs = "scale")%>%
#   ungroup()

#save data to use in fitting the models




dataList <- list(ageAtHarvestData = ageAtHarvestData,
                 sprawningData = sprawningAtHarvestData,
                 indData = indDataWithCatchments,
                 #popnData = popnData,
                 lengthAtAgeData = lengthAtAgeData
                 )

save(dataList, 
     file = "dataset/formattedDataList.RData")

# check the catchment variables and correlation with other variables

#lengthAtAgeThisYear
#GGally::ggpairs(ageAtHarvestData[, c(1:3, 15, 36:45)])

# Catchment variables
#GGally::ggpairs(ageAtHarvestData[, c(53:68)])


