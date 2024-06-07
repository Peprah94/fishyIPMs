# Load libraries
library(ggplot2)
library(dplyr)
library(latex2exp)
# Set the theme of the ggplot
theme_set(theme_bw())

# Load data
load("result/fittedModelWithOverdispersion.RData")

# Load formatted data
load("dataset/formattedDataList.RData")

results <- fishModelMCMCrun$summary$all.chains

# Variable selection probabilities
varSelectionProb <- results[ rownames(results)[grepl("psi", rownames(results))], ]%>%
  data.frame()%>%
  mutate(covariate = c("year", "forest", "pastures",
                       "popnDensity", "summerTemp", "summerSnowDepth",
                       "forestsq", "pasturessq", "popnDensitysq",
                       "summerTempsq", "summerSnowDepthsq", "WPUE"),
         posNeg = ifelse(Mean < 0.4 | Mean > 0.6, "postive", "negative")
  )


ggplot(varSelectionProb, aes(x = covariate, y = Mean, colour = posNeg))+
  geom_point()+
  geom_errorbar(aes(ymin = X95.CI_low, 
                    ymax = X95.CI_upp),
                width = .2)+
  scale_color_manual("posNeg", breaks=c("postive", "negative"),values=c("#D55E00", "black"))+
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "red")+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.6))+
  coord_flip() +
  theme(legend.position="none")+
  ylab("Variable Selection Probablity")
  

# Covariate effects
covEffects <- results[ rownames(results)[grepl("bsurv", rownames(results))], ]%>%
  data.frame()

covEffects[2:13,]%>%
  mutate(covariate = c( "year", "forest", "pastures",
                       "popnDensity", "summerTemp", "summerSnowDepth",
                       "forestsq", "pasturessq", "popnDensitysq",
                       "summerTempsq", "summerSnowDepthsq", "WPUE"),
         posNeg = ifelse(X95.CI_low < 0 & X95.CI_upp > 0, "postive", "negative")
  )%>%
  ggplot(., aes(x = covariate, y = Mean, colour = posNeg))+
  geom_point()+
  geom_errorbar(aes(ymin = Mean - St.Dev., 
                    ymax = Mean + St.Dev.),
                width = .2)+
  scale_color_manual("posNeg", breaks=c("postive", "negative"),values=c("black", "#D55E00"))+
  #geom_hline(yintercept = 0, linetype = "dashed", col = "red")+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.6))+
  coord_flip()+
  theme(legend.position="none")+
  ylab("Effect")


# Effects of age, sex and lake on survival probability
randomEffects <- results[ rownames(results)[grepl("sd", rownames(results))], ]%>%
  data.frame()%>%
  mutate(covariate = c("Age", "Lake", "Sex"))%>%
  ggplot(., aes(x = covariate, y = Mean))+
  geom_point()+
  geom_errorbar(aes(ymin = Mean - St.Dev., 
                    ymax = Mean + St.Dev.),
                width = .2)+
  coord_flip()


# Lake effects
lakeEffect <- results[ rownames(results)[grepl("bsurvLake", rownames(results))], ]%>%
  data.frame()

lakeEffect <- lakeEffect[-nrow(lakeEffect), ]%>%
  mutate(lakeID = levels(as.factor(dataList$ageAtHarvestData$InnsjoNr)),
         posNeg = ifelse(X95.CI_low < 0 & X95.CI_upp > 0, "postive", "negative"))%>%
  ggplot(., aes(x = lakeID, y = Mean, colour = posNeg))+
  geom_point()+
  geom_errorbar(aes(ymin = X95.CI_low, 
                    ymax = X95.CI_upp),
                width = .2)+
  scale_color_manual("posNeg", breaks=c("postive", "negative"),values=c("black", "#D55E00"))+
  coord_flip() +
  theme(legend.position="none")+
  ylab(TeX(sprintf("Mean $\\pm 95$ CI")))
  #xlab(expression(Mean %+-% 95 CI))

# Population growth rate by sex and year
popnGrowthRate <- results[ rownames(results)[grepl("lambda", rownames(results))], ] %>%
  cbind(dataList$ageAtHarvestData[,1:3], .) %>%
  data.frame()%>%
  dplyr::group_by(year, sex)%>%
  summarise_at(., c("Mean", "St.Dev.", "X95.CI_low", "X95.CI_upp"), "mean")%>%
  ggplot(., aes(x = year, y = Mean, color = as.factor(sex)))+
  geom_point()+
  geom_line()

# Population growth rate  and year
popnGrowthRate <- results[ rownames(results)[grepl("lambda", rownames(results))], ] %>%
  cbind(dataList$ageAtHarvestData[,1:3], .) %>%
  data.frame()%>%
  dplyr::group_by(year)%>%
  summarise_at(., c("Mean", "St.Dev.", "X95.CI_low", "X95.CI_upp"), "mean")%>%
  ggplot(., aes(x = as.factor(year), y = Mean))+
  geom_point()+
  geom_line()

# Spawning probability
results[ rownames(results)[grepl("spawn", rownames(results))], ]
