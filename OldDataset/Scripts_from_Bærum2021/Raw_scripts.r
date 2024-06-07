# Trondheimdata --> Lifshistorie analyse
Fisk<-read.csv("OrretRoyeAlleMedGarndyp_20110304.csv", sep=";")
orret<-read.csv("~/Orretdata_phd/Finstad_data/orret.csv", sep=";")
Fisk2<-read.csv("Catch_ved_alder_sum.csv", sep=";")
Fisk_hunn<-read.csv("Catch_ved_alder_sum_hunner.csv", sep=";")
library(lattice)
library(FSA)
library(lme4)
library(nlme)
library(mgcv)
library(bbmle)
library(popbio)
library(MuMIn)
library(plyr)
library(ggplot2)
library(AICcmodavg)
##########################################################
#            Lengdefordeling etc                         #

# Sjekker lengdefordelingen pr vann
orret$alder<-as.factor(orret$alder)
orret.hunn<-subset(orret, subset=kjonn==2)
test.est<-NULL
for(i in 1:max(orret.hunn$loknr))
{
tempo<-subset(orret.hunn, subset=loknr==21)
test.lengde.alder<-lm(log(tempo$Lengde)~tempo$alder)
test.lengde.alder$coef
}

#visualisering
for(i in levels(orret.hunn$loknr))
{
tempo<-subset(orret.hunn, subset=loknr==i)
plot(log(tempo$Lengde)~tempo$alder)
}




##########################################################
#                 Vekst-analyse (growth)                 #
##########################################################
vekst<-read.csv("tilvekstdata_temp_nedbor3.csv", sep=";")
#sjekker fil og respons
#vekst<-subset(vekst,vekstalder !="0")
hist(vekst$inst_vekst)

#Sjekker lengde ved alderm for de forskjellige populasjonene

g0 <- ggplot(data1, aes(x=vekstalder, y=LengdeValder, group= klima, colour=klima)) + stat_smooth(span = 0.9, method="lm") + geom_point()+theme_bw()

#div plott
plot(CPUE~tilvekst, data=vekst )
plot(tilvekst~avg_temp, data=vekst )
plot(tilvekst~sum_nedbor, data=vekst )

#vekstmodeller for alle ?r -->lme4
names(vekst)
vekst$kjonn<-as.factor(vekst$kjonn)
vekst$ID<-as.factor(vekst$ID)
vekst$klima<-as.factor(vekst$klima)
vekst$vekstaar<-as.factor(vekst$vekstaar)
vekst$loknr<-as.factor(vekst$loknr)

#subsample of fish age 1 and up
vekst_1_og_opp<-vekst[!vekst$vekstalder==0,] #Remove zero years old
vekst_1_og_opp<-vekst_1_og_opp[!vekst_1_og_opp$alder_aar==vekst_1_og_opp$alder,] #remove last registrated year growth
vekst_1_og_opp<-vekst_1_og_opp[vekst_1_og_opp$LengdeValder>50,] #Remove som strange insidents where 1 year old fish is calculated to be less than 50 mm
vekst_1_og_opp<-vekst_1_og_opp[!vekst_1_og_opp$LengdeValder==10,]


#Subsample age 1
vekst_1<-vekst[vekst$vekstalder==0,]

#Standarise length for each age-class
vekst_1_og_opp<-transform(vekst_1_og_opp, LengdeVvekst_std=ave(LengdeVvekst, vekstalder, FUN=scale))
#Setter NaN til 0, dvs fikser en bug fra scale funksjonen
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

vekst_1_og_opp[is.nan(vekst_1_og_opp)] <- 0
vekst_1_og_opp<-vekst_1_og_opp[complete.cases(vekst_1_og_opp$kjonn), ]

#Adding column with 1. year growth to "vekst_1_og_opp"

vekst_1_og_opp$first_growth<-vekst_1$LengdeValder[match(vekst_1_og_opp$ID, vekst_1$ID)]

##############
# Model with lmer and glmmtmb
##############
library(bbmle)
library(glmmTMB)

#model selection
select=c("nedbor_mai","nedbor_juni","nedbor_juli")
vekst_1$nebor_summer<-rowSums(vekst_1[select])
# Size at age 1
mod1.1.lmer<-lmer(LengdeValder ~ 1 + (1 |loknr),data=vekst_1) #dispformula=~0 forces variance into the random effects
mod1.2.lmer<-lmer(LengdeValder ~ tempratur + (1 |loknr),data=vekst_1) #dispformula=~0 forces variance into the random effects
mod1.3.lmer<-lmer(LengdeValder ~ NaoV + (1 |loknr),data=vekst_1) #dispformula=~0 forces variance into the random effects
mod1.4.lmer<-lmer(LengdeValder ~ nebor_summer + (1 |loknr),data=vekst_1) #dispformula=~0 forces variance into the random effects
mod1.5.lmer<-lmer(LengdeValder ~ nebor_summer  + NaoV + (1 |loknr),data=vekst_1)
mod1.6.lmer<-lmer(LengdeValder ~ nebor_summer  + tempratur + (1 |loknr),data=vekst_1)
mod1.7.lmer<-lmer(LengdeValder ~ NaoV  + tempratur + (1 |loknr),data=vekst_1)
mod1.8.lmer<-lmer(LengdeValder ~ nebor_summer + tempratur + NaoV + (1 |loknr),data=vekst_1)
mod1.9.lmer<-lmer(LengdeValder ~ tempratur * nebor_summer   + (1 |loknr),data=vekst_1)
mod1.10.lmer<-lmer(LengdeValder ~ tempratur* nebor_summer + NaoV + (1 |loknr),data=vekst_1)
mod1.11.lmer<-lmer(LengdeValder ~ NaoV  * tempratur + (1 |loknr),data=vekst_1)
mod1.12.lmer<-lmer(LengdeValder ~ NaoV  * tempratur + nebor_summer + (1 |loknr),data=vekst_1)

AICtab(mod1.1.lmer,mod1.2.lmer,mod1.3.lmer,mod1.4.lmer,mod1.5.lmer,mod1.6.lmer,mod1.7.lmer,mod1.8.lmer,mod1.9.lmer,mod1.10.lmer,mod1.11.lmer,mod1.12.lmer)
plot(fitted(mod1.11.lmer),resid(mod1.11.lmer))

plot(vekst$poly(alder_aar,2),vekst$LengdeValder)

#R2
r.squaredGLMM(mod1.11.lmer)
anova(mod1.7.lmer)
summary(mod1.3.lmer)




# Size at age 1+ and up

#AICctab(mod.1TMB,mod.2TMB,mod.3TMB,mod.4TMB,mod.5TMB,mod.7TMB)
select=c("nedbor_mai","nedbor_juni","nedbor_juli")
vekst_1_og_opp$nebor_summer<-rowSums(vekst_1_og_opp[select])

data1<-vekst_1_og_opp[complete.cases(vekst_1_og_opp$first_growth), ]

mod.1.lmer<-lmer(LengdeValder ~ alder_aar + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
mod.2.lmer<-lmer(LengdeValder ~ alder_aar+CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
#mod.3.lmer<-lmer(LengdeValder ~ alder_aar*CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
AICtab(mod.1.lmer,mod.2.lmer)
mod.4.lmer<-lmer(LengdeValder ~ first_growth*alder_aar+CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
mod.5.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
#mod.5.2.lmer<-lmer(LengdeValder ~ first_growth*alder_aar+CPUE*alder_aar + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
#mod.5.3.lmer<-lmer(LengdeValder ~ first_growth*alder_aar+CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects
#mod.5.4.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + (1 |loknr/ID),data=data1) #dispformula=~0 forces variance into the random effects

AICtab(mod.1.lmer,mod.2.lmer,mod.3.lmer,mod.4.lmer,mod.5.lmer)

mod.6.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + tempratur + (1 |loknr/ID),data=data1)
mod.7.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + nebor_summer + (1 |loknr/ID),data=data1)
mod.8.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + NaoV + (1 |loknr/ID),data=data1)
mod.9.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + nebor_summer  + NaoV + (1 |loknr/ID),data=data1)
mod.10.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + nebor_summer  + tempratur + (1 |loknr/ID),data=data1)
mod.11.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + NaoV  + tempratur + (1 |loknr/ID),data=data1)
mod.12.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + nebor_summer + NaoV + tempratur + (1 |loknr/ID),data=data1)
mod.13.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + tempratur * nebor_summer + (1 |loknr/ID),data=data1)
mod.14.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + tempratur* nebor_summer + NaoV + (1 |loknr/ID),data=data1)
mod.15.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + NaoV  * tempratur + (1 |loknr/ID),data=data1)
mod.16.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + NaoV  * tempratur + nebor_summer + (1 |loknr/ID),data=data1)
AICtab(mod.1.lmer,mod.2.lmer,mod.4.lmer,mod.5.lmer,mod.6.lmer,mod.7.lmer,mod.8.lmer,mod.9.lmer,mod.10.lmer,mod.11.lmer,mod.12.lmer,mod.13.lmer,mod.14.lmer,mod.15.lmer,mod.16.lmer)

mod.17.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + NaoV  * tempratur + alder_aar*tempratur+(1 |loknr/ID),data=data1)
mod.18.lmer<-lmer(LengdeValder ~ first_growth+alder_aar+CPUE + alder_aar*tempratur + (1 |loknr/ID),data=data1)
mod.19.lmer<-lmer(LengdeValder ~ alder_aar+CPUE + alder_aar*tempratur + NaoV  * tempratur + (1 |loknr/ID),data=data1)
mod.20.lmer<-lmer(LengdeValder ~ alder_aar*first_growth+alder_aar+CPUE + NaoV  * tempratur +alder_aar*tempratur + (1 |loknr/ID),data=data1)

mod.21.lmer<-lmer(LengdeValder ~ first_growth+CPUE + NaoV  * tempratur+NaoV  * alder_aar+ alder_aar*tempratur+(1 |loknr/ID),data=data1,na.action = "na.fail")



library(MuMIn)
#dredged.mod.lmer <- dredge(mod.21.lmer)

mod.list.delta2<-list(mod.21.lmer,mod.17.lmer)

mod.21.lmer<-model.avg(mod.list.delta2, fit=T)


AICtab(mod.1.lmer,mod.2.lmer,mod.4.lmer,mod.5.lmer,mod.6.lmer,mod.7.lmer,mod.8.lmer,mod.9.lmer,mod.10.lmer,mod.11.lmer,mod.12.lmer,mod.13.lmer,mod.14.lmer,mod.15.lmer,mod.16.lmer,mod.17.lmer,mod.18.lmer,mod.20.lmer,mod.19.lmer,mod.21.lmer)
require(kimisc)
model_list <- tibble::lst(mod.1.lmer,mod.2.lmer,mod.4.lmer,mod.5.lmer,mod.6.lmer,mod.7.lmer,mod.8.lmer,mod.9.lmer,mod.10.lmer,mod.11.lmer,mod.12.lmer,mod.13.lmer,mod.14.lmer,mod.15.lmer,mod.16.lmer,mod.17.lmer,mod.18.lmer,mod.20.lmer,mod.19.lmer,mod.21.lmer)


# Compare models with AIC table
aic_table <- aictab(model_list)

T# Add a column to the table with variable names
model_vars <- sapply(model_list, function(x) gsub(":", "*", paste(names(x$means), collapse = " + ")))
aic_table <- cbind(Modnames = names(model_vars), model_vars) %>%
  merge(aic_table, by="Modnames") %>%
  arrange(AICc)

AICtab(mod.1.lmer,mod.2.lmer,mod.3.lmer,mod.4.lmer,mod.5.lmer,mod.5.2.lmer,mod.6.lmer,mod.7.lmer,mod.8.lmer,mod.9.lmer,mod.10.lmer,mod.11.lmer,mod.12.lmer,mod.13.lmer,mod.14.lmer,mod.15.lmer,mod.16.lmer,mod.17.lmer,mod.18.lmer,mod.20.lmer,mod.19.lmer,mod.21.lmer,mod.22.lmer,mod.23.lmer,mod.24.lmer,mod.25.lmer,mod.26.lmer)



#sjekk
plot(fitted(mod.21.lmer),resid(mod.21.lmer))


#R2
r.squaredGLMM(mod.21.lmer)
anova(mod.17.lmer)
summary(mod.21.lmer)

#plot(mod.2TMB)

## Prediksjoner
# Str alder 1
#mod1.13.lmer<-lmer(LengdeValder ~ CPUE + tempratur * NaoV + (1 |loknr),data=vekst_1)
pred_frame1 <- with(vekst_1, expand.grid(tempratur=seq(min(tempratur),max(tempratur), by=1),NaoV= seq(-1.5,1.5, by=0.2), loknr=levels(as.factor(loknr))))

pred<-predict(mod1.11.lmer,pred_frame1,
              se.fit = F)

pred_frame1$pred<-as.data.frame(pred)[,1]


first_year_growth<-ggplot(pred_frame1, aes(tempratur, NaoV)) + geom_raster(aes(fill = pred)) +
  #scale_fill_gradientn(name = "Swim speed",colours=c("white","yellow","red"),limits=c(0, 0.68), breaks=seq(0,0.68,by=0.1))+
  scale_fill_viridis_c(name = "First year growth (mm)",limits=c(43, 85), breaks=seq(45,85,by=5),option = "C")+scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))+
  guides(fill = guide_colorbar(barheight = 1.5,barwidth = 45,title.vjust =0.25,label.position = "top",frame.colour = "black", ticks.colour = "black"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14))+
  ylab(expression ("Winter NAO"))+
  xlab(expression("Summer temperature"))+
  theme(axis.text=element_text(size=18,colour="black"),axis.title=element_text(size=18,face="bold")) +
  theme(axis.text.y = element_text(angle = 90,hjust=0.2),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(legend.position="top")


#contourplot(pred~tempratur * NaoV, data=pred_frame1, label.style="align", xlab="Temp", ylab="NaoV",cuts=30, as.table=T )
#trellis.focus("panel", 1, 1)
#lpoints(vekst_1$tempratur,vekst_1$NaoV)
#trellis.unfocus()

# Str fra alder 1 og opp
#mod.1.lmer<-lmer(LengdeValder ~ LengdeVvekst_std*vekstalder + CPUE  + NaoV + (1 |loknr/ID),data=vekst_1_og_opp) #dispformula=~0 forces variance into the random effects

pred_frame1 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=seq(min(NaoV),max(NaoV),by=0.2),tempratur=max(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame2 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=seq(min(NaoV),max(NaoV),by=0.2),tempratur=min(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame3 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=max(NaoV),tempratur=seq(min(tempratur),max(tempratur),by=1),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame4 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=min(NaoV),tempratur=seq(min(tempratur),max(tempratur),by=1),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame5 <- with(data1, expand.grid(first_growth=45,CPUE=mean(CPUE),NaoV=seq(min(NaoV),max(NaoV),by=0.2),tempratur=seq(min(tempratur),max(tempratur),by=1),alder_aar=seq(min(alder_aar),10, by=0.5),ID=1,loknr=levels(factor(loknr))))
pred_frame6 <- with(data1, expand.grid(first_growth=85,CPUE=mean(CPUE),NaoV=seq(min(NaoV),max(NaoV),by=0.2),tempratur=seq(min(tempratur),max(tempratur),by=1),alder_aar=seq(min(alder_aar),10, by=0.5),ID=1,loknr=levels(factor(loknr))))

pred_frame7 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=max(NaoV),tempratur=max(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame8 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=min(NaoV),tempratur=max(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame9 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=max(NaoV),tempratur=min(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))
pred_frame10 <- with(data1, expand.grid(first_growth=seq(min(45),max(85), by=2),CPUE=mean(CPUE),NaoV=min(NaoV),tempratur=min(tempratur),alder_aar=seq(2,10, by=1),ID=1,loknr=levels(factor(loknr))))

pred<-predict(mod.21.lmer,as.data.frame(pred_frame1),allow.new.levels=T,
              se.fit = F)
pred2<-predict(mod.21.lmer,as.data.frame(pred_frame2),allow.new.levels=T,
              se.fit = F)
pred3<-predict(mod.21.lmer,as.data.frame(pred_frame3),allow.new.levels=T,
              se.fit = F)
pred4<-predict(mod.21.lmer,as.data.frame(pred_frame4),allow.new.levels=T,
               se.fit = F)
pred5<-predict(mod.21.lmer,as.data.frame(pred_frame5),allow.new.levels=T,
               se.fit = F)
pred6<-predict(mod.21.lmer,as.data.frame(pred_frame6),allow.new.levels=T,
               se.fit = F)

pred7<-predict(mod.21.lmer,as.data.frame(pred_frame7),allow.new.levels=T,
               se.fit = F)
pred8<-predict(mod.21.lmer,as.data.frame(pred_frame8),allow.new.levels=T,
               se.fit = F)
pred9<-predict(mod.21.lmer,as.data.frame(pred_frame9),allow.new.levels=T,
               se.fit = F)
pred10<-predict(mod.21.lmer,as.data.frame(pred_frame10),allow.new.levels=T,
               se.fit = F)

#pred3<-predict(mod.30.lmer,as.data.frame(pred_frame3),allow.new.levels=T,
             # se.fit = F)
#pred4<-predict(mod.30.lmer,as.data.frame(pred_frame4),allow.new.levels=T,
              # se.fit = F)


pred_frame1$pred<-as.data.frame(pred)[,1]
pred_frame2$pred<-as.data.frame(pred2)[,1]
pred_frame3$pred<-as.data.frame(pred3)[,1]
pred_frame4$pred<-as.data.frame(pred4)[,1]
pred_frame5$pred<-as.data.frame(pred5)[,1]
pred_frame6$pred<-as.data.frame(pred6)[,1]

pred_frame7$pred<-as.data.frame(pred7)[,1]
pred_frame8$pred<-as.data.frame(pred8)[,1]
pred_frame9$pred<-as.data.frame(pred9)[,1]
pred_frame10$pred<-as.data.frame(pred10)[,1]

varying_nao <- ggplot(pred_frame3, aes(x=alder_aar, y=pred)) +stat_smooth(colour="black",size=1,level=0.99)+theme_bw()+stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame4,colour="black",size=1,level=0.99,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("Growth (mm)")+
  xlab("Age")+
  labs(title = "Winter NAO")+
  theme(plot.title = element_text(size=18,hjust=1))+
  coord_cartesian(xlim=c(1, 10))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

varying_temp <- ggplot(pred_frame1, aes(x=alder_aar, y=pred)) +stat_smooth(colour="black",size=1,level=0.99)+theme_bw()+stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame2,colour="black",size=1,level=0.99,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("Length (mm)")+
  xlab("Age")+
  labs(title = "Summer temperature")+
  theme(plot.title = element_text(size=18,hjust=1))+
  coord_cartesian(xlim=c(1, 10))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))


varying_first_year <- ggplot(pred_frame5, aes(x=alder_aar, y=pred)) +stat_smooth(colour="black",size=1,level=0.99)+theme_bw()+stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame6,colour="black",size=1,level=0.99,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("Length (mm)")+
  xlab("Age")+
  labs(title = "First year growth")+
  theme(plot.title = element_text(size=18,hjust=1))+
  coord_cartesian(xlim=c(1, 10))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))


varying_temp_and_nao <- ggplot(pred_frame7, aes(x=alder_aar, y=pred)) +stat_smooth(colour="black",size=1,level=0.99)+theme_bw()+
  stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame8,colour="black",size=1,level=0.99,linetype="dotdash")+
  stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame9,colour="red",size=1,level=0.99,linetype="dotdash")+
 stat_smooth(aes(x=alder_aar, y=pred),data=pred_frame10,colour="blue",size=1,level=0.99,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("Length (mm)")+
  xlab("Age")+
  labs(title = "NAO and summer temperature")+
  theme(plot.title = element_text(size=14,hjust=1))+
  coord_cartesian(xlim=c(1, 10, by=1))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

pred_frame_growth<-cbind(pred_frame7,pred_frame8$pred,pred_frame9$pred,pred_frame10$pred)

require(reshape2)
require(data.table)
pred_frame_growth_long <- melt(setDT(pred_frame_growth), id.vars = c("first_growth","CPUE","NaoV","tempratur","alder_aar"),
                               measure.vars =  c("pred","pred_frame8$pred","pred_frame9$pred","pred_frame10$pred"),
                        variable.name = 'Temperature_scenario', value.name = "Length")
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

tgc <- summarySE(pred_frame_growth_long, measurevar="Length", groupvars=c("alder_aar","Temperature_scenario"))





varying_temp_and_nao <- ggplot(tgc, aes(x=alder_aar, y=Length, colour=Temperature_scenario,fill=Temperature_scenario))+theme_bw()+
  geom_ribbon(aes(x=,alder_aar, ymax = Length+sd,ymin=Length-sd),
                                colour=NA, alpha = 0.2)+
  stat_smooth(method =lm,size=1.5,show.legend=FALSE)+
  stat_smooth(method =lm,fill=NA)+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14))+
  scale_fill_discrete(name="Climate scenario", labels = c("High w-NAO, high s-temp", "Low w-NAO, high s-temp", "High w-NAO, low s-temp","Low w-NAO, low s-temp"))+
  scale_colour_discrete(name="Climate scenario", labels = c("High w-NAO, high s-temp", "Low w-NAO, high s-temp", "High w-NAO, low s-temp","Low w-NAO, low s-temp"))+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("Length (mm)")+
  #labs(title = "NAO and summer temperature")+
  theme(plot.title = element_text(size=18,hjust=1))+
  scale_x_continuous(name="Age", breaks =seq(2, 10,1))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

library(stargazer)
stargazer(mod.21.lmer,mod1.11.lmer,
          type="html",
          title="Regression Results",
          intercept.bottom = F,
          intercept.top = T,
          ci = F, digits=2,
          notes = "This is a caption.",
          model.names = T,
          single.row = T,
          out="table_growth.doc")

##########################################################
#                    Modningsanalyse                     #
##########################################################
(orret$moden)
#maturation ogives
names(orret)

#orret.hunn$kjonn<-as.factor(orret$kjonn)
orret.hunn$loknr<-as.factor(orret.hunn$loknr)
orret.hunn$moden1<-as.numeric(orret.hunn$moden1)

names(orret)
table(orret.hunn$alder_aar) # Trimmer bort alder 10 11 og 12, tilsammen 13 individer
orret.hunn<-(subset(orret.hunn, subset=alder_aar!=12))
orret.hunn<-(subset(orret.hunn, subset=alder_aar!=11))
#orret<-(subset(orret, subset=alder_aar!=10))
# NB! Dette gjorde s? bergna R2 sank fra 40% til 23%!
plot(orret$moden1, orret$lengde)
# Setter loknr som random-effect. M? sjekke dette ut mer.
#orret<-subset(Fisk, art==1)
orret.alder<-subset(orret.hunn, alder_aar!="NA")
orret.alder$loknr<-as.factor(orret.alder$loknr)
orret.alder$kjonn<-as.factor(orret.alder$kjonn)

#tar ut populasjoner med ingen eller f? (dvs < 10) modne individer
keep.pop<-NULL
for(i in 1:21)
{
  temp.pop<-subset(orret.alder, subset=loknr==i)
  if(dim(subset(temp.pop, subset=moden1==1))[1]<5) next
  keep.pop<-rbind(keep.pop,temp.pop)
}



#inkluderer moh i fullfil
####################nao<-read.csv("nao.csv",header=T, sep=";")
#keep.pop$moh<-NULL
#for(i in 1:length(keep.pop$klima)){
#  tmp=subset(nao,Loknr==keep.pop$loknr[i])
#  keep.pop$moh[i]=tmp$moh
#}

#keep.popF<-subset(keep.pop, kjonn=="2")

#includerer snitt-temp for de forskjellige lokalitetene
keep.pop$mean_temp<-1
keep.pop$loknr<-factor(keep.pop$loknr)
for(i in levels(keep.pop$loknr)){
  tmp=vekst_1_og_opp[vekst_1_og_opp$loknr==i,]
  keep.pop[keep.pop$loknr==i,]$mean_temp=mean(tmp$tempratur)
}

# includes variation
keep.pop$var_temp<-1
for(i in levels(keep.pop$loknr)){
  tmp=vekst_1_og_opp[vekst_1_og_opp$loknr==i,]
  keep.pop[keep.pop$loknr==i,]$var_temp=var(tmp$tempratur)
}


#includerer precipitation for de forskjellig
keep.pop$nebor_summer<-1
for(i in levels(keep.pop$loknr)){
  tmp=vekst_1[vekst_1$loknr==i,]
  keep.pop[keep.pop$loknr==i,]$nebor_summer=mean(tmp$nedbor)
}

#includerer cpue for de forskjellig
keep.pop$CPUE<-1
for(i in levels(keep.pop$loknr)){
  tmp=vekst_1[vekst_1$loknr==i,]
  keep.pop[keep.pop$loknr==i,]$CPUE=mean(tmp$CPUE)
}


#lager modeller
prob.spawn1<-glmer(moden1~alder_aar+log(Lengde) + (1|loknr), data=keep.pop, family=binomial)
prob.spawn2<-glmer(moden1~alder_aar*log(Lengde) + (1|loknr), data=keep.pop, family=binomial)

AICctab(prob.spawn1,prob.spawn2)

keep.pop$var_temp_std <- (keep.pop$var_temp - mean(keep.pop$var_temp)) / sd(keep.pop$var_temp)

hist(keep.pop$var_temp_std)
#Effect of var_temp
#prob.spawn3<-glmer(moden1~alder_aar*log(Lengde)+var_temp + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn4<-glmer(moden1~alder_aar*log(Lengde) + log(Lengde)*var_temp+ (1|loknr), data=keep.pop, family=binomial)
#prob.spawn5<-glmer(moden1~alder_aar*var_temp+log(Lengde)*var_temp + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn6<-glmer(moden1~alder_aar*log(Lengde) +alder_aar*var_temp + (1|loknr), data=keep.pop, family=binomial)


#Effect of mean_temp
prob.spawn9<-glmer(moden1~alder_aar*log(Lengde)+log(mean_temp) + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn8<-glmer(moden1~alder_aar*log(Lengde) + log(Lengde)*log(mean_temp) +(1|loknr), data=keep.pop, family=binomial)
#prob.spawn5<-glmer(moden1~alder_aar*var_temp+log(Lengde)*var_temp + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn9<-glmer(moden1~alder_aar*log(Lengde) +alder_aar*log(mean_temp) + (1|loknr), data=keep.pop, family=binomial)




#Effect of nedbor_summer
prob.spawn4<-glmer(moden1~alder_aar*log(Lengde)+log(nebor_summer) + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn11<-glmer(moden1~alder_aar*log(Lengde) + log(Lengde)*log(nebor_summer) +(1|loknr), data=keep.pop, family=binomial)
#prob.spawn5<-glmer(moden1~alder_aar*var_temp+log(Lengde)*var_temp + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn12<-glmer(moden1~alder_aar*log(Lengde) +alder_aar*log(nebor_summer) + (1|loknr), data=keep.pop, family=binomial)

#effect of density
prob.spawn5<-glmer(moden1~alder_aar*log(Lengde)+CPUE + (1|loknr), data=keep.pop, family=binomial)
#prob.spawn6<-glmer(moden1~alder_aar*log(Lengde)+log(mean_temp) + var_temp  + (1|loknr), data=keep.pop, family=binomial)

AICctab(prob.spawn1,prob.spawn2,prob.spawn9,prob.spawn4,prob.spawn5)


#summary(prob.spawnT)
#visualisering LAG PLOTT SOM IKKE ER SOMA INN; DVS MED SAMME SKALA. LEGG OGS? P? POLYNOMER
#orret2<-subset(orret, alder_aar==0:13)
     data_kald<-with(keep.pop,expand.grid(moden1=0,alder_aar=seq(1, 10, by=1),Lengde=seq(min(Lengde), max(Lengde), by=10),mean_temp=min(mean_temp),loknr=0))
     data_varm<-with(keep.pop,expand.grid(moden1=0,alder_aar=seq(1, 10, by=1),Lengde=seq(min(Lengde), max(Lengde), by=10),mean_temp=max(mean_temp),loknr=0))
     #data2<-with(b,expand.grid(moden1=0,alder_aar=seq(1, 10, by=1),Lengde=seq(min(Lengde), max(Lengde), by=10),moh=seq(0, 100,by=10),LengdeValder=seq(45,65,by=1), kjonn=2,loknr=0))

     data_kald$pred.mature <- plogis(predict(prob.spawn9, data_kald,allow.new.levels =T))
     data_varm$pred.mature <- plogis(predict(prob.spawn9, data_varm,allow.new.levels =T))
     data_varm[data_varm$alder_aar==1,]$pred.mature<-0
     data_varm[data_varm$Lengde<150,]$pred.mature<-0
     grob_kald <- grobTree(textGrob("Cold summers", x=1,  y=420, hjust=0,
                                    gp=gpar(col="black", fontsize=14, fontface="bold")))

     grob_varm <- grobTree(textGrob("Warm summers", x=0.1,  y=450, hjust=0,
                                    gp=gpar(col="black", fontsize=14, fontface="bold")))

     tgc$alder_aar<-as.numeric(tgc$alder_aar)
     names(tgc)[names(tgc) == "Length"] <- "Lengde"
     modning2011<-ggplot(data_kald, aes(alder_aar, Lengde)) + geom_raster(aes(fill = pred.mature)) +
       scale_fill_gradientn(name = "Probability for spawning",colours=c("white","green"),limits=c(0, 1), breaks=seq(0,1,by=0.1))+
       scale_x_continuous(expand = c(0, 0),breaks =round(seq(1,20,by=1),1)) +
       scale_y_continuous(expand=c(0, 0),breaks =round(seq(100,450,by=50),1))+
       guides(fill = guide_colorbar(barheight = 1.5,barwidth = 45,title.vjust =0.25,label.position = "top",frame.colour = "black", ticks.colour = "black"))+
       theme(legend.title = element_text(size=18, face="bold"))+
       theme(legend.text = element_text(size=18))+
       ylab("Length (mm)")+
       xlab("Age")+
       theme(axis.text=element_text(size=24,colour="black"),axis.title=element_text(size=24,face="bold")) +
       theme(axis.text.y = element_text(angle = 90,hjust=0.2),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
       theme(legend.position="top")+
       geom_smooth(data=tgc, aes(x=alder_aar, y=Lengde),method="lm",formula = y ~ splines::bs(x, 3),size=1.5)+
       geom_segment(aes(x = 6, y = 255, xend = 0.6, yend = 255), color="black",linetype="dotted",size=1)  +
       geom_segment(aes(x = 6, y = 90, xend = 6, yend = 255), color="black",linetype="dotted",size=1)+
       geom_text(x=5, y=400, label="Low s-temp",size=8,colour="black")



     modning2017<-ggplot(data_varm, aes(alder_aar, Lengde)) + geom_raster(aes(fill = pred.mature)) +
       scale_fill_gradientn(name = "Probability for spawning",colours=c("white","green"),limits=c(0, 1), breaks=seq(0,1,by=0.1))+
       theme(axis.text.y.right = element_text(angle = 270, hjust = 0.5))+
       scale_x_continuous(expand = c(0, 0),breaks =round(seq(1,20,by=1),1)) +
       scale_y_continuous(expand=c(0, 0),breaks =round(seq(100,450,by=50),1),position = "right")+
       guides(fill = guide_colorbar(barheight = 1.5,barwidth = 45,title.vjust =0.25,label.position = "top",frame.colour = "black", ticks.colour = "black"))+
       theme(legend.title = element_text(size=18, face="bold"))+
       theme(legend.text = element_text(size=18))+
       ylab("Length (mm)")+
       xlab("Age")+
       theme(axis.text=element_text(size=24,colour="black"),axis.title=element_text(size=24,face="bold")) +
       theme(axis.text.y = element_text(angle = 90,hjust=0.2),panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
       theme(legend.position="top")+
       geom_smooth(data=tgc, aes(x=alder_aar, y=Lengde),method="lm",formula = y ~ splines::bs(x, 3),size=1.5)+
       geom_segment(aes(x = 4.55, y = 200, xend = 10.4, yend = 200), color="black",linetype="dotted",size=1)  +
       geom_segment(aes(x = 4.55, y = 90, xend = 4.55, yend = 200), color="black",linetype="dotted",size=1)+
       geom_text(x=5, y=400, label="High s-temp",size=8,colour="black")

     ggpubr::ggarrange(modning2011, modning2017, ncol=2, nrow=1, common.legend = TRUE, legend="top")








##########################################################
#                 d?delighetsanalyse                     #
##########################################################
Catch_hunn<-subset(orret.alder, kjonn==2)
orret.alder$loknr<-as.factor(orret.alder$loknr)
Catch_hunn$klima<-as.factor(Catch_hunn$klima)

#Making dataset where slopes are decreasing for each lake
diff_set <- function(x, gap) {
  ind <- c(F, diff(x) > gap)
  if(sum(ind) == 0) return(!ind)
  subst <- x[-unique(c(which(ind), which(ind)-1))]
  x %in% subst
}
catch.curve.data=NULL
#for(i in levels(orret.alder$loknr)){
#tmp=subset(orret.alder,loknr==i)
catch.curve.data=ddply(orret.alder,~alder_aar,summarise,sum=sum(s_vekting))
#a$loknr=rep(as.numeric(tmp$loknr[1]),length(levels(as.factor(a$alder_aar))))
#a<-as.data.frame(a)
catch.curve.data<-catch.curve.data[catch.curve.data$alder_aar>4,]
#b<-b[diff_set(b$sum, 0),]
#catch.curve.data<-rbind(b,catch.curve.data)
#}






chapmanRobson(sum~alder_aar,data=catch.curve.data)

#$est
#Estimate Std. Error
#S 48.8188553 1.29591064
#Z  0.7163522 0.05161701
#keep.Z<-cbind(keep.Z,cr$est)
#}
#keep.Z

#A=1???e???Z

survival<-(2.71828182846^-0.7008791)

# Visualisering Z og S
Fisk_hunn2$loknr<-as.factor(Fisk_hunn2$loknr)
boxplot(Z~loknr, data=Fisk_hunn)
names(Fisk_hunn)
   # ...med ggplot2
ggplot(Fisk_hunn, aes(x=loknr, y=Z)) +
geom_errorbar(aes(ymin=Z-Std_err_Z, ymax=Z+Std_err_Z), width=.5) +
geom_point()


Fisk_hunn$lengde<-as.numeric(Fisk_hunn$lengde)
Fisk_hunn$loknr<-as.factor(Fisk_hunn$loknr)
#
lengdetest<-lm(Fisk_hunn$alder~log(Fisk_hunn$lengde)*Fisk_hunn$loknr)
test.data<-subset(orret3, subset=loknr==10)
pooled.mort.data<-data.frame(age=seq(1,8,1), catch=as.numeric(tapply(test.data$s_vekting,as.factor(test.data$alder),sum)))
cc.p<-with(pooled.mort.data,catchCurve(age,catch,5:10))
coef(cc.p)
confint(cc.p)
cc<-with(subset(mort.data.adj, subset=year==2009),catchCurve(age,catch,4:10))
coef(cc)
confint(cc)


           mort.mod<-lm(log(catch)~age*tot.catch, data=subset(mort.data.adj, subset=age>3))
           growth.season<-data.frame(year=c(2006,2007,2008,2009),gr.L=c(122,90,105,114))
           plot(growth.season$gr.L*as.numeric(tapply(mort.data.adj$ma2.tot.catch,as.factor(mort.data.adj$year),mean)),keep.S)
           pred.mort<-growth.season$gr.L*as.numeric(tapply(mort.data.adj$ma2.tot.catch,as.factor(mort.data.adj$year),mean))
           fit.S<-lm(as.numeric(keep.S)~pred.mort)
           ann.S<-data.frame(year=c(2006,2007,2008,2009),S=as.numeric(keep.S))


           mort.mod<-lm(LnN~Egentlig_alder*N, data=subset(Fisk_hunn, subset=Egentlig_alder>3))


# Modelerer S med klima #
#########################

#Datasjekk
par(mfrow=c(2,2))
plot(Fisk_hunn$S)
boxplot(S, data=Fisk_hunn)
hist(Fisk_hunn$S, main="")
qqnorm(Fisk_hunn$S)

names(Fisk_hunn)
boxplot(Fisk_hunn$S~Fisk_hunn$klima)

Fisk_hunn$klima<-as.factor(Fisk_hunn$klima)
test.S<-lmer(S~klima+gjen_L_4aar + (1|loknr), data=Fisk_hunn)
test.S2<-lmer(S~klima+gjen_L_4aar + (0 + gjen_L_4aar|loknr), data=Fisk_hunn)

#glm-tiln?rming
glm.Z<-glm(lnN~alder, data=Fisk_hunn)
glm.Z.1<-glm(LnN~poly(alder,2)*as.factor(klima), data=Fisk_hunn)
glm.Z.2<-glm(LnN~alder*gjen_L_4aar, data=Fisk_hunn)
glm.Z.3<-glm(LnN~alder+gjen_L_4aar, data=Fisk_hunn)



summary(glm.Z.3)

summary(as.factor(orret$alder))


ggplot(Fisk_hunn, aes(x=Alder, y=lnN)) +stat_smooth(colour="black",size=1,level=0.99)+theme_bw()+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  ylab("lnN")+
  xlab("Age")+
  labs(title = "C")+
  theme(plot.title = element_text(size=18,hjust=1))+
  coord_cartesian(xlim=c(4.5, 9),ylim =c(1,5) )+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

           ###################
           #Populations matix#
           ###################


#########################################################################
#                   Genererer 100-?rsklimasenario          #
#########################################################################
##UTEN TEMPTREND, sampler tilfeldig innenfor temp de siste 10 ?r)

lambda.vecutrend=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:100){
temp.utrend<-as.vector(rnorm(100,11.78,1.32))
NAO.utrend<-as.vector(runif(100,-1.5,1.5))
n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue, fordel p? 4 til siste levende
CPUE.sim<-c(3.06) #CPUE (snitt)
L.start<-c(50,92,129,159,184,185,186,187,188,189)
keep.lengde<-NULL
keep.lengde.vecutrend<-NULL
lambda.vec=NULL
pop.matrix=NULL

for(i in 1:100) {


  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
 L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,mean_temp=mean(temp.utrend),loknr=1)
  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
  sex.ratio<-1
    s1 <- srFuns("Shepherd")
    n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
    #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
    surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
    exp.eggs<-n.eggs*surv1
    surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

}

lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)
temp.utrend=NULL
keep.lengde.vecutrend<-c(keep.lengde,keep.lengde.vecutrend)
setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend1<-cbind(lambda.vecutrend,simu,year,n.init)
lambda.vecutrend1<-data.frame(lambda.vecutrend1)



####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
lambda.vecmtrend=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:100){

temp.100=NULL
simu=NULL
year=NULL
temp.utrend<-NULL
temp.mtrend<-NULL
temp.utrend<-as.vector(rnorm(100,11.78,1.32))
NAO.utrend<-as.vector(runif(100,-1.5,1.5))

for (i in 1:100){
  s<-temp.utrend[i]
  e<-s+(i*0.04)
  temp.mtrend<-c(temp.mtrend,e)
}



n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
CPUE.sim<-c(3.06) #CPUE (snit
L.start<-c(50,92,129,159,184,185,186,187,188,189)
keep.lengde<-NULL
keep.lengde.vecmtrend<-NULL
tilv.i=NULL
L.a=NULL
L.a.i=NULL

lambda.vec=NULL
pop.matrix=NULL

for(i in 1:100) {


  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,mean_temp=mean(temp.utrend)+(i*0.04),loknr=1)
  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
  sex.ratio<-1
  s1 <- srFuns("Shepherd")
  n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
  #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
  surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
  exp.eggs<-n.eggs*surv1
  surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

}

lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
temp.mtrend=NULL
setTxtProgressBar(pb, g)
keep.lengde.vecmtrend<-c(keep.lengde,keep.lengde.vecmtrend)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend1<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend1<-data.frame(lambda.vecmtrend1)

#############################################
#temp uten trend men mer og mer upredikabelt
lambda.vecupred=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:100){
temp.100=NULL
simu=NULL
year=NULL
temp.utrend<-NULL
temp.utrend.upred<-NULL
temp.utrend<-as.vector(rnorm(100,11.78,1.32))
#NAO.utrend<-as.vector(runif(100,-1.5,1.5))
for (i in 1:100){
  s<-temp.utrend[i]
  e<-s+(i*rnorm(1, mean = 0, sd = 0.03))
  temp.utrend.upred<-c(temp.utrend.upred,e)
}

NAO.utrend<-as.vector(runif(100,-1.5,1.5))

n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
CPUE.sim<-c(3.06) #CPUE (snitt)
L.start<-c(50,92,129,159,184,185,186,187,188,189)
keep.lengde<-NULL
keep.lengde.vecupred<-NULL
tilv.i=NULL
L.a=NULL
L.a.i=NULL
lambda.vec=NULL
pop.matrix=NULL

for(i in 1:100) {


  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,mean_temp=mean(temp.utrend.upred),loknr=1)
  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
  sex.ratio<-1
  s1 <- srFuns("Shepherd")
  n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
  #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
  surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
  exp.eggs<-n.eggs*surv1
  surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

}
#n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
#CPUE.sim<-c(3.06) #CPUE (snitt)
#L.start<-c(50,92,129,159,184,185,186,187,188,189)
#keep.lengde<-NULL
#tilv.i=NULL
#L.a=NULL
#L.a.i=NULL

#lambda.vec=NULL
#pop.matrix=NULL

#for(i in 1:100) {


#  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
#  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
#  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
#  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
#  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
#  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
#  a$index<-rep(i,10)
#  keep.lengde<-rbind(keep.lengde,a)
#  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


 # L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
#  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
#  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
#  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
#  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
#  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
#  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
#  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
#  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
#  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
#  tilv.i=NULL
#  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

#  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
#  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
#  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
#  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
#  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
#  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
#  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
#  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
#  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
#  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
#  L.start<-L.a.i



#  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a,mean_temp=mean(temp.utrend),loknr=1)
#  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
#  sex.ratio<-1
#  s1 <- srFuns("Shepherd")
#  n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
  #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
#  surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
#  exp.eggs<-n.eggs*surv1
#  surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

#  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
#  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
#  proj.t<- pop.projection(A.trout, n.init, 2)


#  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
#  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
#  n.init<-proj.t$stage.vectors[,2]
#  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

#}
lambda.vecupred<-c(lambda.vec,lambda.vecupred)
setTxtProgressBar(pb, g)
temp.utrend.upred=NULL
keep.lengde.vecupred<-c(keep.lengde,keep.lengde.vecupred)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupred)/100))
year<-rep(1:100, 100)
lambda.vecupred1<-cbind(lambda.vecupred,simu,year)
lambda.vecupred1<-data.frame(lambda.vecupred1)





#Med temptrend men mer og mer upredikabelt
temp.100=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
lambda.vecupredmtrend=NULL
for (g in 1:100){
temp.utrend<-NULL
temp.utrend<-as.vector(rnorm(100,11.78,1.32))
temp.mtrend.upred<-NULL
    for (i in 1:100){
  s<-temp.utrend[i]
  b<-s+(i*0.02)
  e<-b+(i*rnorm(1, mean = 0, sd = 0.03))
  temp.mtrend.upred<-c(temp.mtrend.upred,e)
}


NAO.utrend<-as.vector(runif(100,-1.5,1.5))

n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
CPUE.sim<-c(3.06) #CPUE (snitt)
L.start<-c(50,92,129,159,184,185,186,187,188,189)
keep.lengde<-NULL
keep.lengde.vecupredmtrend<-NULL
tilv.i=NULL
L.a=NULL
L.a.i=NULL
lambda.vec=NULL
pop.matrix=NULL

for(i in 1:100) {


  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,mean_temp=mean(temp.utrend)+(i*0.02),loknr=1)
  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
  sex.ratio<-1
  s1 <- srFuns("Shepherd")
  n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
  #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
  surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
  exp.eggs<-n.eggs*surv1
  surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

}
lambda.vecupredmtrend<-c(lambda.vec,lambda.vecupredmtrend)
temp.mtrend.upred=NULL
setTxtProgressBar(pb, g)
keep.lengde.vecupredmtrend<-c(keep.lengde,keep.lengde.vecupredmtrend)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupredmtrend)/100))
year<-rep(1:100, 100)
lambda.vecupredmtrend1<-cbind(lambda.vecupredmtrend,simu,year)
lambda.vecupredmtrend1<-data.frame(lambda.vecupredmtrend1)





#Visualisering med ggplot
g1 <- ggplot(lambda.vecutrend1, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw()+
geom_smooth(aes(x=year, y=lambda.vecupredmtrend),data=lambda.vecupredmtrend1,linetype="dotdash",colour="red",size=1)+
geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend1,size=1,colour="red")+
  geom_smooth(aes(x=year, y=lambda.vecupred),data=lambda.vecupred1,colour="blue",linetype="dotdash",size=1)+
  labs(title="Historical NAO trend") + theme(plot.title = element_text(vjust=-7, hjust = 0.9))+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(0, 100),ylim=c(0.985,1.07))+
  ylab("Lambda ")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

#Boxplot
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}



pred_frame<-cbind(lambda.vecutrend1,lambda.vecupredmtrend1$lambda.vecupredmtrend,lambda.vecmtrend1$lambda.vecmtrend,lambda.vecupred1$lambda.vecupred)
pred_frame<-pred_frame[pred_frame$year>50,]
require(reshape2)
require(data.table)
pred_frame_long <- melt(setDT(pred_frame), measure.vars =  c("lambda.vecutrend", "lambda.vecmtrend1$lambda.vecmtrend","lambda.vecupred1$lambda.vecupred","lambda.vecupredmtrend1$lambda.vecupredmtrend"),
                        variable.name = 'Temperature_scenario', value.name = "Lambda")

tgc <- summarySE(pred_frame_long, measurevar="Lambda", groupvars=c("Temperature_scenario"))

ggplot(tgc, aes(x=Temperature_scenario, y=Lambda,group=Temperature_scenario,colour=Temperature_scenario,)) +theme_bw()+
  geom_linerange(aes(ymin=Lambda-sd, ymax=Lambda+sd),size=1.5) +
  #geom_line(aes(colour=Perception_of_statment), size=1.5) +
  #scale_x_discrete(breaks=c("Uenig","Verken_eller","Enig"),
  #                 labels=c("Disagree", "Not sure", "Agree"),expand = c(0.05,0.05))+
  geom_point(size=2) +
  #scale_colour_discrete(name="Perception of statement", labels = c("Dont know", "Guessing", "Manipulation","Political","Research"))+
  #facet_grid( .~gender)+
  ylab("Predicted probability")+
  xlab("Trust in carnivore science")+
  #labs(title = "A")+
  #theme_minimal()+
  theme(legend.title = element_text(size=18, face="bold"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  #theme(legend.title = element_text(size=18, face="bold"),legend.position="top")+
  theme(legend.text = element_text(size=18))+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(colour="lightgrey"),
    panel.grid.minor = element_blank(),
    #  ,panel.border = element_blank()
    legend.background = element_blank(),
    #legend.key = element_rect(fill = "white", color = NA),
    # Change legend key size and key width
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(3.5,"cm"),
    #axis.title.y  = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
  )



#########################################################################
#                   Genererer 100-?rsklimasenario med h?y NAO            #
#########################################################################
##UTEN TEMPTREND, sampler tilfeldig innenfor temp de siste 10 ?r)

lambda.vecutrend=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:100){
  temp.utrend<-as.vector(rnorm(100,11.78,1.32))
  NAO.utrend<-as.vector(rnorm(100,0.5,0.5))
  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue, fordel p? 4 til siste levende
  CPUE.sim<-c(3.06) #CPUE (snitt)
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL

  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {


    pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
    pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
    a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a,mean_temp=mean(temp.utrend),loknr=1)
    prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


    #fecunditet
    sex.ratio<-1
    s1 <- srFuns("Shepherd")
    n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
    #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
    surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
    exp.eggs<-n.eggs*surv1
    surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



    #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }

  lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)
  temp.utrend=NULL
  #keep.lengde.vecutrend<-c(keep.lengde,keep.lengde.vecutrend)
  setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend2<-cbind(lambda.vecutrend,simu,year)
lambda.vecutrend2<-data.frame(lambda.vecutrend2)



####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
lambda.vecmtrend=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:100){

  temp.100=NULL
  simu=NULL
  year=NULL
  temp.utrend<-NULL
  temp.mtrend<-NULL
  temp.utrend<-as.vector(rnorm(100,11.78,1.32))
  NAO.utrend<-as.vector(rnorm(100,0.5,0.5))

  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*0.04)
    temp.mtrend<-c(temp.mtrend,e)
  }



  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
  CPUE.sim<-c(3.06) #CPUE (snit
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL
  tilv.i=NULL
  L.a=NULL
  L.a.i=NULL

  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {


    pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
    pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
    a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a,mean_temp=mean(temp.utrend)+(i*0.04),loknr=1)
    prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


    #fecunditet
    sex.ratio<-1
    s1 <- srFuns("Shepherd")
    n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
    #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
    surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
    exp.eggs<-n.eggs*surv1
    surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



    #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }

  lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
  temp.mtrend=NULL
  setTxtProgressBar(pb, g)
  #keep.lengde.vecmtrend<-c(keep.lengde,keep.lengde.vecmtrend)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend2<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend2<-data.frame(lambda.vecmtrend2)


#temp uten trend men mer og mer upredikabelt
lambda.vecupred=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:10){
  temp.100=NULL
  simu=NULL
  year=NULL
  temp.utrend<-NULL
  temp.utrend.upred<-NULL
  temp.utrend<-as.vector(rnorm(100,11.78,1.32))
  #NAO.utrend<-as.vector(runif(100,-1.5,1.5))
  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*rnorm(1, mean = 0, sd = 0.03))
    temp.utrend.upred<-c(temp.utrend.upred,e)
  }

  NAO.utrend<-as.vector(rnorm(100,0.5,0.5))

  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
  CPUE.sim<-c(3.06) #CPUE (snitt)
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL
  keep.lengde.vecupred<-NULL
  tilv.i=NULL
  L.a=NULL
  L.a.i=NULL
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {


    pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
    pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
    a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,mean_temp=mean(temp.utrend.upred),loknr=1)
    prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


    #fecunditet
    sex.ratio<-1
    s1 <- srFuns("Shepherd")
    n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
    #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
    surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
    exp.eggs<-n.eggs*surv1
    surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



    #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }
  #n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
  #CPUE.sim<-c(3.06) #CPUE (snitt)
  #L.start<-c(50,92,129,159,184,185,186,187,188,189)
  #keep.lengde<-NULL
  #tilv.i=NULL
  #L.a=NULL
  #L.a.i=NULL

  #lambda.vec=NULL
  #pop.matrix=NULL

  #for(i in 1:100) {


  #  pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.utrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
  #  pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
  #  pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
  #  pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
  #  a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
  #  a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
  #  a$index<-rep(i,10)
  #  keep.lengde<-rbind(keep.lengde,a)
  #  pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


  # L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  #  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  #  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  #  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  #  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  #  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  #  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  #  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  #  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  #  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  #  tilv.i=NULL
  #  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  #  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  #  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  #  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  #  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  #  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  #  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  #  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  #  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  #  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  #  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  #  L.start<-L.a.i



  #  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a,mean_temp=mean(temp.utrend),loknr=1)
  #  prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


  #fecunditet
  #  sex.ratio<-1
  #  s1 <- srFuns("Shepherd")
  #  n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
  #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
  #  surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
  #  exp.eggs<-n.eggs*surv1
  #  surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



  #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

  #  A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
  #  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  #  proj.t<- pop.projection(A.trout, n.init, 2)


  #  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  #  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  #  n.init<-proj.t$stage.vectors[,2]
  #  CPUE.sim<-sum(n.init[4:10])

  # update progress bar

  #}
  lambda.vecupred<-c(lambda.vec,lambda.vecupred)
  setTxtProgressBar(pb, g)
  temp.utrend.upred=NULL
  keep.lengde.vecupred<-c(keep.lengde,keep.lengde.vecupred)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupred)/100))
year<-rep(1:100, 100)
lambda.vecupred2<-cbind(lambda.vecupred,simu,year)
lambda.vecupred2<-data.frame(lambda.vecupred2)






#Med temptrend men mer og mer upredikabelt
temp.100=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
lambda.vecupredmtrend=NULL
for (g in 1:100){
  temp.utrend<-NULL
  temp.utrend<-as.vector(rnorm(100,11.78,1.32))
  temp.mtrend.upred<-NULL
  for (i in 1:100){
    s<-temp.utrend[i]
    b<-s+(i*0.02)
    e<-b+(i*rnorm(1, mean = 0, sd = 0.03))
    temp.mtrend.upred<-c(temp.mtrend.upred,e)
  }


  NAO.utrend<-as.vector(rnorm(100,0.5,0.5))

  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
  CPUE.sim<-c(3.06) #CPUE (snitt)
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL
  tilv.i=NULL
  L.a=NULL
  L.a.i=NULL
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {


    pframe1aar<-with(vekst_1_og_opp,expand.grid(LengdeValder=0,first_growth=0,alder_aar=seq(1,10,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend.upred[i],NaoV=NAO.utrend[i], ID=1,loknr=1, vekstaar=0))
    pframe1aar$first_growth<-predict(mod1.11.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$Lengde <- predict(mod.21.lmer,pframe1aar,allow.new.levels =T)
    pframe1aar$alder_aar<-as.factor(pframe1aar$alder_aar)
    a=ddply(pframe1aar,~alder_aar,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a[a$alder_aar==1,]$mean<-mean(pframe1aar$first_growth)
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$alder_aar<-as.numeric(pframe1aar$alder_aar)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a,mean_temp=mean(temp.utrend)+(i*0.02),loknr=1)
    prob.mat<-plogis(predict(prob.spawn9,data.t,allow.new.levels =T))


    #fecunditet
    sex.ratio<-1
    s1 <- srFuns("Shepherd")
    n.eggs <-exp(log(data.t$Lengde*2.21-6.15))*prob.mat*sex.ratio # how many offspring permature adult (asumes mature if prob is more than 0.5)
    #n.eggs <- round(n.eggs*0.5) # Asumes 50% females, and just reduces the number of eggs by 50%
    surv1 <- exp(log(s1(sum(n.eggs,na.rm = T),a=0.04,b=0.000003,c=3.5))-log(sum(n.eggs,na.rm = T)))
    exp.eggs<-n.eggs*surv1
    surv<-c(survival*0.60,survival,survival*1.2,survival*1.2,survival*1.2,survival,survival,survival,survival,0)



    #exp.reprod<-n.eggs#fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.eggs[3],exp.eggs[4],exp.eggs[5],exp.eggs[6],exp.eggs[7],exp.eggs[8],exp.eggs[9],exp.eggs[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }
  lambda.vecupredmtrend<-c(lambda.vec,lambda.vecupredmtrend)
  temp.mtrend.upred=NULL
  setTxtProgressBar(pb, g)
  #keep.lengde.vecupredmtrend<-c(keep.lengde,keep.lengde.vecupredmtrend)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupredmtrend)/100))
year<-rep(1:100, 100)
lambda.vecupredmtrend2<-cbind(lambda.vecupredmtrend,simu,year)
lambda.vecupredmtrend2<-data.frame(lambda.vecupredmtrend2)


#Visualisering med ggplot
g2 <- ggplot(lambda.vecutrend2, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw()+
  geom_smooth(aes(x=year, y=lambda.vecupredmtrend),data=lambda.vecupredmtrend2,linetype="dotdash",colour="red",size=1)+
  geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend2,size=1,colour="red")+
  geom_smooth(aes(x=year, y=lambda.vecupred),data=lambda.vecupred2,colour="blue",linetype="dotdash",size=1)+
  labs(title="Increased NAO trend") + theme(plot.title = element_text(vjust=-7, hjust = 0.9))+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(0, 100),ylim=c(0.985,1.07))+
  ylab("Lambda")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))

figure <- ggarrange(g1, g2,
                    ncol = 1, nrow = 2)


pred_frame2<-cbind(lambda.vecutrend2,lambda.vecupredmtrend2$lambda.vecupredmtrend,lambda.vecmtrend2$lambda.vecmtrend,lambda.vecupred2$lambda.vecupred)
pred_frame2<-pred_frame2[pred_frame2$year>50,]
pred_frame_long2 <- melt(setDT(pred_frame2), measure.vars =  c("lambda.vecutrend", "lambda.vecmtrend2$lambda.vecmtrend","lambda.vecupred2$lambda.vecupred","lambda.vecupredmtrend2$lambda.vecupredmtrend"),
                        variable.name = 'Temperature_scenario', value.name = "Lambda")

tgc2 <- summarySE(pred_frame_long2, measurevar="Lambda", groupvars=c("Temperature_scenario"))



#########################################################################
#                   Genererer 100-?rsklimasenario for sone 3            #
#########################################################################

temp100<-subset(vekst, klima==3)
pb<-txtProgressBar(min = 0, max = 1000, style = 3)

lambda.vecutrend=NULL
for (g in 1:1000){
temp.utrend=NULL
temp.utrend<-as.vector(runif(100,9.4,12.56))
CPUE.sim=3.15
n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
L.start<-c(45,85,126,165,201,235,267,295,345,350)
keep.lengde<-NULL


surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0.61)
lambda.vec=NULL
pop.matrix=NULL
for(i in 1:100) {
  pframe1aar<-with(vekst,expand.grid(LengdeValder=0,alder_aar=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$alder_aar==0,0,NA)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$alder_aar==1,L.start[1],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$alder_aar==2,L.start[2],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$alder_aar==3,L.start[3],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
  mm <- model.matrix(terms(lme.lengde),pframe1aar)
  pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
  pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
  a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep(800,10))
  mm.t <- model.matrix(terms(prob.spawnT),data.t)
  prob.mat<- mm.t %*% fixef(prob.spawnT)


  #fecunditet

  fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
  keep.fec<-fec

  sex.ratio<-0.5
  exp.reprod<-fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

}

lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)

setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend3<-cbind(lambda.vecutrend,simu,year)
lambda.vecutrend3<-data.frame(lambda.vecutrend3)



####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
lambda.vecmtrend<-NULL
for (g in 1:1000){
temp.100=NULL
simu=NULL
year=NULL
temp.mtrend=NULL
temp.utrend<-NULL
temp.utrend<-as.vector(runif(100,9.4,12.56))

  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*0.04)
    temp.mtrend<-c(temp.mtrend,e)
  }

CPUE.sim=3.15
n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
L.start<-c(45,85,126,165,201,235,267,295,345,350)
keep.lengde<-NULL


surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0.61)
lambda.vec=NULL
pop.matrix=NULL
for(i in 1:100) {
  pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
  mm <- model.matrix(terms(lme.lengde),pframe1aar)
  pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
  pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
  a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i




  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep(800,10))
  mm.t <- model.matrix(terms(prob.spawnT),data.t)
  prob.mat<- mm.t %*% fixef(prob.spawnT)


  #fecunditet

  fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
  keep.fec<-fec

  sex.ratio<-0.5
  exp.reprod<-fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

}

lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
setTxtProgressBar(pb, g)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend3<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend3<-data.frame(lambda.vecmtrend3)


#temp uten trend men mer og mer upredikabelt
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
lambda.vecupred<-NULL
for (g in 1:1000){
temp.100=NULL
temp.utrend.upred=NULL
simu=NULL
year=NULL
temp.utrend<-NULL
temp.utrend<-as.vector(runif(100,9.4,12.56))

  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*rnorm(1, mean = 0, sd = 0.05))
    temp.utrend.upred<-c(temp.utrend.upred,e)
  }






CPUE.sim=3.15
n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
L.start<-c(45,85,126,165,201,235,267,295,345,350)
keep.lengde<-NULL


surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0.61)
lambda.vec=NULL
pop.matrix=NULL
for(i in 1:100) {
  pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.utrend.upred[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
  mm <- model.matrix(terms(lme.lengde),pframe1aar)
  pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
  pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
  a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i



  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,kjonn=2,moh=rep(800,10))
  mm.t <- model.matrix(terms(prob.spawnT),data.t)
  prob.mat<- mm.t %*% fixef(prob.spawnT)


  #fecunditet

  fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
  keep.fec<-fec

  sex.ratio<-0.5
  exp.reprod<-fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])

}
lambda.vecupred<-c(lambda.vec,lambda.vecupred)
setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupred)/100))
year<-rep(1:100, 100)
lambda.vecupred3<-cbind(lambda.vecupred,simu,year)
lambda.vecupred3<-data.frame(lambda.vecupred3)





#Med temptrend men mer og mer upredikabelt
lambda.vecupredmtrend=NULL
pb<-txtProgressBar(min = 0, max = 1000, style = 3)
for (g in 1:1000){
temp.100=NULL
simu=NULL
year=NULL
temp.utrend<-NULL
temp.utrend<-as.vector(runif(100,9.4,12.56))
temp.mtrend.upred<-NULL

  for (i in 1:100){
    s<-temp.utrend[i]
    b<-s+(i*0.04)
    e<-b+(i*rnorm(1, mean = 0, sd = 0.05))
    temp.mtrend.upred<-c(temp.mtrend.upred,e)
  }




CPUE.sim=3.15
n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
L.start<-c(45,85,126,165,201,235,267,295,345,350)
keep.lengde<-NULL


surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0)
lambda.vec=NULL
pop.matrix=NULL
for(i in 1:100){
  pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend.upred[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
  pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
  mm <- model.matrix(terms(lme.lengde),pframe1aar)
  pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
  pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
  a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
  a$index<-rep(i,10)
  keep.lengde<-rbind(keep.lengde,a)
  pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


  L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
  L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
  L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
  L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
  L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
  L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
  L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
  L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
  L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
  L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
  tilv.i=NULL
  tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

  L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
  L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
  L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
  L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
  L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
  L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
  L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
  L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
  L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
  L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
  L.start<-L.a.i


  data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,kjonn=2,moh=rep(800,10))
  mm.t <- model.matrix(terms(prob.spawnT),data.t)
  prob.mat<- mm.t %*% fixef(prob.spawnT)


  #fecunditet

  fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
  keep.fec<-fec

  sex.ratio<-0.5
  exp.reprod<-fec*sex.ratio*plogis(prob.mat)

  A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
  diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
  proj.t<- pop.projection(A.trout, n.init, 2)


  lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
  pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
  n.init<-proj.t$stage.vectors[,2]
  CPUE.sim<-sum(n.init[4:10])
}
lambda.vecupredmtrend<-c(lambda.vec,lambda.vecupredmtrend)
setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecupredmtrend)/100))
year<-rep(1:100, 100)
lambda.vecupredmtrend3<-cbind(lambda.vecupredmtrend,simu,year)
lambda.vecupredmtrend3<-data.frame(lambda.vecupredmtrend3)




#Visualisering med ggplot

  g3 <- ggplot(lambda.vecutrend3, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw() +
  geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend3,colour="black",size=1,linetype="dotdash")+
  geom_smooth(aes(x=year, y=lambda.vecupred),data=lambda.vecupred3,colour="black",size=1)+
  geom_smooth(aes(x=year, y=lambda.vecupredmtrend),data=lambda.vecupredmtrend3,colour="red",size=1)+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(45, 100),ylim=c(0.99,1.005))+
  ylab("Lambda ")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))


par(ylog=TRUE, cex=1.2)
plot(1:100,seq(1,100,by=100), type="n", xlab="Year", ylab="Number of individuals",log="y")
for(i in 1:10)
{
  lines(1:100,pop.matrix[i,], col=i)
}





par(ylog=TRUE, cex=1.2)
plot(1:100,seq(1,20000, by=200), type="n", xlab="Year", ylab="Number of individuals",log="y")
for(i in 1:10)
{
  lines(1:100,pop.matrix[i,], col=i)
}
text(80,10,paste("Lambda =",round(mean(lambda.vec),3),"?",round(sd(lambda.vec),3)))





##############################
# Scenarior hvor modning g?r mot lavereliggende moh
##############################

#########################################################################
#                   Genererer 100-?rsklimasenario for sone 1            #
#########################################################################
##UTEN TEMPTREND, sampler tilfeldig innenfor temp de siste 10 ?r)


temp100<-subset(vekst, klima==1) #legg til nedb?r
lambda.vecutrend=NULL
pb<-txtProgressBar(min = 0, max = 100, style = 3)
for (g in 1:100){
  temp.utrend<-as.vector(runif(100,12,15))
  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue, fordel p? 4 til siste levende
  CPUE.sim<-c(3.06) #CPUE (snitt)
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL

  surv<-c(0.034,0.70,0.65,0.47,0.37,0.37,0.37,0,0,0)
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],nedbor=seq(358, 509, by=30),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(10, 150,by=10), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep(50,10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }

  lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)
  temp.utrend=NULL
  setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend1<-cbind(lambda.vecutrend,simu,year)
lambda.vecutrend1<-data.frame(lambda.vecutrend1)



####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
lambda.vecmtrend=NULL
pb<-txtProgressBar(min = 0, max = 100, style = 3)
for (g in 1:100){

  temp.100=NULL
  simu=NULL
  year=NULL
  temp.utrend<-NULL
  temp.mtrend<-NULL
  temp.utrend<-as.vector(runif(1000,12,15))


  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*0.04)
    temp.mtrend<-c(temp.mtrend,e)
  }



  n.init<-c(21.47,2.14,1.72,1.29,0.86,0.43,0.21,0.15,0.086,0.032)#lag fordeling i forhold til cpue
  CPUE.sim<-c(3.06) #CPUE (snit
  L.start<-c(50,92,129,159,184,185,186,187,188,189)
  keep.lengde<-NULL
  tilv.i=NULL
  L.a=NULL
  L.a.i=NULL
  surv<-c(0.034,0.70,0.65,0.47,0.37,0.37,0.37,0,0,0)
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],nedbor=seq(358, 509, by=30),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(10, 150,by=10), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep((50-i),10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])
    # update progress bar

  }

  lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
  temp.mtrend=NULL
  setTxtProgressBar(pb, g)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend1<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend1<-data.frame(lambda.vecmtrend1)






#Visualisering med ggplot
g1 <- ggplot(lambda.vecutrend1, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw()+
  #geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend1,colour="black",size=1,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(45, 100),ylim=c(0.8,1.2))+
  ylab("Lambda ")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))





#########################################################################
#                   Genererer 100-?rsklimasenario for sone 2            #
#########################################################################
##UTEN TEMPTREND, sampler tilfeldig innenfor temp de siste 10 ?r)


temp100<-subset(vekst, klima==2)
lambda.vecutrend=NULL
pb<-txtProgressBar(min = 0, max = 100, style = 3)
for (g in 1:100){

  temp.100=NULL
  temp.utrend<-NULL
  temp.utrend<-as.vector(runif(100,10.6,13.52))
  CPUE.sim=7.2
  n.init<-c(50.5,5.0,4.04,3.03, 2.02,1.01,0.5,0.35,0.2,0.076)
  L.start<-c(41,80,119,155,197,242,255,265,295,325)
  keep.lengde<-NULL


  surv<-c(0.035,0.70,0.65,0.46,0.46,0.46,0.46,0.46,0.46,0.46)
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],nedbor=seq(383, 491, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(160, 400,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i





    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep(300,10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

    # update progress bar

  }

  lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)
  temp.utrend=NULL
  setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend2<-cbind(lambda.vecutrend,simu,year)
lambda.vecutrend2<-data.frame(lambda.vecutrend2)





####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
lambda.vecmtrend=NULL
pb<-txtProgressBar(min = 0, max = 100, style = 3)
for(g in 1:100){

  temp.100=NULL
  simu=NULL
  year=NULL
  temp.utrend<-NULL
  temp.mtrend<-NULL
  temp.utrend<-as.vector(runif(100,10.6,13.52))


  for(i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*0.04)
    temp.mtrend<-c(temp.mtrend,e)
  }


  CPUE.sim=7.2
  n.init<-c(50.5,5.0,4.04,3.03, 2.02,1.01,0.5,0.35,0.2,0.076)
  L.start<-c(41,80,119,155,197,242,255,265,295,325)
  keep.lengde<-NULL


  surv<-c(0.035,0.70,0.65,0.46,0.46,0.46,0.46,0.46,0.46,0.46)
  lambda.vec=NULL
  pop.matrix=NULL

  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],nedbor=seq(383, 491, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(160, 400,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i





    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep((300-(i*2)),10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])
  }

  lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
  temp.mtrend=NULL
  setTxtProgressBar(pb, g)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend2<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend2<-data.frame(lambda.vecmtrend2)



#Visualisering med ggplot
g2 <- ggplot(lambda.vecutrend2, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw()+
  geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend2,colour="black",size=1,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(45, 100),ylim=c(0.995,1.01))+
  ylab("Lambda")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))





#########################################################################
#                   Genererer 100-?rsklimasenario for sone 3            #
#########################################################################

temp100<-subset(vekst, klima==3)
pb<-txtProgressBar(min = 0, max = 100, style = 3)

lambda.vecutrend=NULL
for (g in 1:100){
  temp.utrend=NULL
  temp.utrend<-as.vector(runif(100,9.4,12.56))
  CPUE.sim=3.15
  n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
  L.start<-c(45,85,126,165,201,235,267,295,345,350)
  keep.lengde<-NULL


  surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0.61)
  lambda.vec=NULL
  pop.matrix=NULL
  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.utrend[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i



    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep(800,10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

  }

  lambda.vecutrend<-c(lambda.vec,lambda.vecutrend)

  setTxtProgressBar(pb, g)
}
simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecutrend)/100))
year<-rep(1:100, 100)
lambda.vecutrend3<-cbind(lambda.vecutrend,simu,year)
lambda.vecutrend3<-data.frame(lambda.vecutrend3)



####MED TEMPTREND, sampler tilfeldig ?r men ?ker temp med 4 grader p? 100 ?r
pb<-txtProgressBar(min = 0, max = 100, style = 3)
lambda.vecmtrend<-NULL
for (g in 1:100){
  temp.100=NULL
  simu=NULL
  year=NULL
  temp.mtrend=NULL
  temp.utrend<-NULL
  temp.utrend<-as.vector(runif(100,9.4,12.56))

  for (i in 1:100){
    s<-temp.utrend[i]
    e<-s+(i*0.04)
    temp.mtrend<-c(temp.mtrend,e)
  }

  CPUE.sim=3.15
  n.init<-c(22.1,2.21,1.77,1.33,0.88,0.44,0.22,0.15,0.09,0.03)
  L.start<-c(45,85,126,165,201,235,267,295,345,350)
  keep.lengde<-NULL


  surv<-c(0.028,0.50,0.65,0.61,0.61,0.61,0.61,0.61,0.61,0.61)
  lambda.vec=NULL
  pop.matrix=NULL
  for(i in 1:100) {
    pframe1aar<-with(vekst,expand.grid(LengdeValder=0,vekstalder=seq(0,9,by=1),CPUE=CPUE.sim,tempratur=temp.mtrend[i],nedbor=seq(291, 405, by=20),NaoV=seq(min(NaoV),max(NaoV), by=1),moh=seq(500, 900,by=30), kjonn=2, ID=0,loknr=0, vekstaar=0))
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==0,0,NA)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==1,L.start[1],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==2,L.start[2],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==3,L.start[3],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==4,L.start[4],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==5,L.start[5],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==6,L.start[6],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==7,L.start[7],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==8,L.start[8],pframe1aar$LengdeVvekst)
    pframe1aar$LengdeVvekst<-ifelse(pframe1aar$vekstalder==9,L.start[9],pframe1aar$LengdeVvekst)
    mm <- model.matrix(terms(lme.lengde),pframe1aar)
    pframe1aar$Lengde <- mm %*% fixef(lme.lengde)
    pframe1aar$vekstalder<-as.factor(pframe1aar$vekstalder)
    a=ddply(pframe1aar,~vekstalder,summarise,mean=mean(Lengde),sd=sd(Lengde))
    a$index<-rep(i,10)
    keep.lengde<-rbind(keep.lengde,a)
    pframe1aar$vekstalder<-as.numeric(pframe1aar$vekstalder)


    L.a<-subset(keep.lengde,index==i,select=mean)[1:10,1]
    L.a[2]<-ifelse(L.a[2]<L.a[1],L.a[1]+5,L.a[2])
    L.a[3]<-ifelse(L.a[3]<L.a[2],L.a[2]+5,L.a[3])
    L.a[4]<-ifelse(L.a[4]<L.a[3],L.a[3]+5,L.a[4])
    L.a[5]<-ifelse(L.a[5]<L.a[4],L.a[4]+5,L.a[5])
    L.a[6]<-ifelse(L.a[6]<L.a[5],L.a[5]+5,L.a[6])
    L.a[7]<-ifelse(L.a[7]<L.a[6],L.a[6]+5,L.a[7])
    L.a[8]<-ifelse(L.a[8]<L.a[7],L.a[7]+5,L.a[8])
    L.a[9]<-ifelse(L.a[9]<L.a[8],L.a[8]+5,L.a[9])
    L.a[10]<-ifelse(L.a[10]<L.a[9],L.a[9]+5,L.a[10])
    tilv.i=NULL
    tilv.i<-c(L.a[2]-L.start[1],L.a[3]-L.start[2],L.a[4]-L.start[3],L.a[5]-L.start[4],L.a[6]-L.start[5],L.a[7]-L.start[6],L.a[8]-L.start[7],L.a[9]-L.start[8],L.a[10]-L.start[9])

    L.a.i<-c(L.a[1],(L.start[-10]+tilv.i))
    L.a.i[2]<-ifelse(L.a.i[2]<L.a.i[1],L.a.i[1]+1,L.a.i[2])
    L.a.i[3]<-ifelse(L.a.i[3]<L.a.i[2],L.a.i[2]+1,L.a.i[3])
    L.a.i[4]<-ifelse(L.a.i[4]<L.a.i[3],L.a.i[3]+1,L.a.i[4])
    L.a.i[5]<-ifelse(L.a.i[5]<L.a.i[4],L.a.i[4]+1,L.a.i[5])
    L.a.i[6]<-ifelse(L.a.i[6]<L.a.i[5],L.a.i[5]+1,L.a.i[6])
    L.a.i[7]<-ifelse(L.a.i[7]<L.a.i[6],L.a.i[6]+1,L.a.i[7])
    L.a.i[8]<-ifelse(L.a.i[8]<L.a.i[7],L.a.i[7]+1,L.a.i[8])
    L.a.i[9]<-ifelse(L.a.i[9]<L.a.i[8],L.a.i[8]+1,L.a.i[9])
    L.a.i[10]<-ifelse(L.a.i[10]<L.a.i[9],L.a.i[9]+1,L.a.i[10])
    L.start<-L.a.i




    data.t<-data.frame(moden1=rep(0,10),alder_aar=seq(1,10,1),Lengde=L.a.i,moh=rep((800-(i*3)),10))
    mm.t <- model.matrix(terms(prob.spawnT),data.t)
    prob.mat<- mm.t %*% fixef(prob.spawnT)


    #fecunditet

    fec<-exp(log(data.t$Lengde)*2.21-6.15)                  #beregn fecunditet
    keep.fec<-fec

    sex.ratio<-0.5
    exp.reprod<-fec*sex.ratio*plogis(prob.mat)

    A.trout<-matrix(c(0, 0, exp.reprod[3],exp.reprod[4],exp.reprod[5],exp.reprod[6],exp.reprod[7],exp.reprod[8],exp.reprod[9],exp.reprod[10],rep(0,90)), nrow = 10, byrow = TRUE)
    diag(A.trout[2:10,1:9])<-surv[1:9];A.trout[dim(A.trout)[1],dim(A.trout)[2]]<-surv[length(surv)]
    proj.t<- pop.projection(A.trout, n.init, 2)


    lambda.vec<-c(lambda.vec,as.numeric(lambda(A.trout)))
    pop.matrix<-cbind(pop.matrix,proj.t$stage.vectors[,2])
    n.init<-proj.t$stage.vectors[,2]
    CPUE.sim<-sum(n.init[4:10])

  }

  lambda.vecmtrend<-c(lambda.vec,lambda.vecmtrend)
  setTxtProgressBar(pb, g)
}

simu=NULL
year=NULL
simu<-as.factor(ceiling(seq_along(lambda.vecmtrend)/100))
year<-rep(1:100, 100)
lambda.vecmtrend3<-cbind(lambda.vecmtrend,simu,year)
lambda.vecmtrend3<-data.frame(lambda.vecmtrend3)



#Visualisering med ggplot

g32 <- ggplot(lambda.vecutrend3, aes(x=year, y=lambda.vecutrend)) +geom_smooth()+theme_bw() +
  geom_smooth(aes(x=year, y=lambda.vecmtrend),data=lambda.vecmtrend3,colour="black",size=1,linetype="dotdash")+
  theme(legend.background = element_rect(fill="white",colour="white"))+
  coord_cartesian(xlim=c(40, 100),ylim=c(0.995,1.01))+
  ylab("Lambda ")+
  xlab("Year")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"))+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +


  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))
