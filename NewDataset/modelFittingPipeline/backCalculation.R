library(RFishBC)
library(FSA)
library(tidyverse)

setwd("C:/GitHub/fishyIPMs/NewDataset")
fish_data <- read.table("dataset/innlandsfisk.txt", header = TRUE, sep = ";", dec = ".")
focal_speciesid <- 5



species_data<-fish_data[fish_data$FK_ArtID==focal_speciesid,]


#remove NAs for Length and Radius 
backcalc_data<-species_data[complete.cases(species_data$Radius), ]
backcalc_data<-backcalc_data[complete.cases(backcalc_data$Lengde), ]
#Check fish data for different lakes
table(backcalc_data$FK_InnsjoNr)
# Some lakes have very few individuals across time(less than 10), remove these lakes as they dont seem representative for the species
# Remove rows with less than min_rows_per_category per category
backcalc_data <- backcalc_data %>%
  group_by(FK_InnsjoNr) %>%
  filter(n() >= 9) %>%
  ungroup()



plot(backcalc_data$Lengde~backcalc_data$Radius)
#Seem to be two types of measurements, just remove the "outliers" for now
backcalc_data<-backcalc_data[backcalc_data$Radius<400,]
#Cheking for strange length Radius correlations among lakes and years
data_check <- backcalc_data %>%
  group_by(FK_InnsjoNr, Aar) %>%
  summarize(R_squared = summary(lm(Lengde ~ Radius))$r.squared,
            Intercept = coef(lm(Lengde ~ Radius))[1])
print(data_check)

#decide on cut-off on 0.8
ok_lakes_years<-data_check[data_check$R_squared>0.8,]
#Remove r_squared of 1, which seems wronge
ok_lakes_years<-ok_lakes_years[!ok_lakes_years$R_squared==1,]

backcalc_data<-backcalc_data[backcalc_data$FK_InnsjoNr %in% ok_lakes_years$FK_InnsjoNr & backcalc_data$Aar %in% ok_lakes_years$Aar ,]
#plot(backcalc_data$Lengde~backcalc_data$Radius)

#setting up data for RFishBC (renaming variables)
backcalc_data <- backcalc_data %>%
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
backcalc_data_FL <- backcalc_data %>%
  group_by(FK_InnsjoNr) %>%
  do(backCalc(., Lengde, BCM = "FRALE",
              inFormat = "wide", outFormat = "long", digits = 0)) %>%
  ungroup()


# Print the resulting dataframe
print(backcalc_data_FL)

###
#Rename variables
backcalc_data_FL <- backcalc_data_FL %>%
  rename(LengdeValder = bclen,
         alder_aar = ann)

#including year of length increment (specific year affecting specific growth)
backcalc_data_FL$vekstaar<-backcalc_data_FL$Aar-((backcalc_data_FL$agecap-backcalc_data_FL$alder_aar+1))
