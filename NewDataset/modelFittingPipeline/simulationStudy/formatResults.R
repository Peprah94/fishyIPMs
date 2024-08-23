library(dplyr)
library(ggplot2)
library(reshape2)

#Load simulated data

load("NewDataset/result/simulatedDatasets.RData")
load("NewDataset/result/fittedModelWithSimData.RData")

# 
covEffects <- lapply(simModelFit, function(x){
  x[rownames(x)[grepl("bsurv",rownames(x))][1:6], "Median"] - c(-0.05, 0.2, 0, 0.5, -0.8, 0)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt(.)%>%
  dplyr::rename(variable = Var2)%>%
  ggplot()+
  geom_boxplot(mapping = aes(x= variable, y = value))


survProb <- lapply(simModelFit, function(x){
  x[rownames(x)[grepl("bsurvProbAge",rownames(x))][1:6], "Median"] - c(0.9, 0.5, 0.2, 0.1, 0.02, 0.01)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt(.)%>%
  dplyr::rename(variable = Var2)%>%
  ggplot()+
  geom_boxplot(mapping = aes(x= variable, y = value))

spawnProb <- lapply(simModelFit, function(x){
  x[rownames(x)[grepl("spawn",rownames(x))], "Median"] - c(0.1)
})%>%
  do.call("rbind", .)%>%
  reshape2::melt(.)%>%
  dplyr::rename(variable = Var2)%>%
  ggplot()+
  geom_boxplot(mapping = aes(x= variable, y = value))
