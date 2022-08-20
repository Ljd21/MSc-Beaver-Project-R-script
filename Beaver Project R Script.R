##LUCAS DE RUIG THESIS R SCRIPT
#Clear Workspace
rm(list=ls())

#install and download necessary packages
install.packages("readr")
install.packages("ggplot2")
install.packages("plyr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("car")
install.packages("MuMIn")
install.packages("ggeffects")
install.packages("Rmisc")
install.packages("lme4")
install.packages("plotly")
install.packages("devtools")
install_github("vqv/ggbiplot")
library(devtools)
library(tidyverse)
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(car)
library(MuMIn)
library(ggeffects)
library(Rmisc)
library(emmeans)
library(ggbiplot)
library(plotly)

#Set Working Directory
setwd("C:/Users/lderu/Desktop/Beaver Project files/Beaver project excel files")

##BIOASSAY ANALYSIS (OAK AND WETTEX)
#Load data
Bioassay_Data <- read.csv(file = "Bioassay Data.csv")

# Format data
Bioassay_Data$Site<-as.factor(Bioassay_Data$Site)
Bioassay_Data$Location<-as.factor(Bioassay_Data$Location)
Bioassay_Data$Bioassay<-as.factor(Bioassay_Data$Bioassay)
Bioassay_Data$Mesh<-as.factor(Bioassay_Data$Mesh)

#Simplify df and pivot
Bioassay_Data<- dplyr::select(Bioassay_Data, Site, Location, Bioassay, Mesh, k_R1, k_R2, k_R3, k_R4) %>%
  dplyr::rename(R1 = k_R1,
                R2 = k_R2,
                R3 = k_R3,
                R4 = k_R4,
                ) %>%
  pivot_longer(
    cols=starts_with("R"),
    names_to = "Replicate",
    values_to = "k",
    values_drop_na = TRUE)

## ENCLOSURE VS. CONTROL
#Create separate dataframes for Oak and Wettex, Enclosure and Control
Bioassay_Data_EC <- dplyr::filter(Bioassay_Data, Location == "Enclosure" | Location == "Control")
Oak_EC <- dplyr::filter(Bioassay_Data_EC, Bioassay == "Oak")
Wettex_EC <- dplyr::filter(Bioassay_Data_EC, Bioassay == "Wettex")  

#LMM Enclosure and Control (EC)
#Oak EC
Lmer_oak_EC <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Oak_EC)
summary(Lmer_oak_EC)
Anova(Lmer_oak_EC, test.statistic = "F")
#check model assumptions
plot(Lmer_oak_EC)
#qqplot
qqnorm(resid(Lmer_oak_EC))
qqline(resid(Lmer_oak_EC))
# transform k values
Oak_EC$Kt <- log10(Oak_EC$k)
# refit LMM
Lmer_oak_EC <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Oak_EC, na.action = na.omit)
summary(Lmer_oak_EC)
Anova(Lmer_oak_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_EC)
#qqplot
qqnorm(resid(Lmer_oak_EC))
qqline(resid(Lmer_oak_EC))
#explained variation by model
r.squaredGLMM(Lmer_oak_EC)

#Oak Fine EC
Oak_f_EC <- filter(Oak_EC, Mesh == "Fine")
Lmer_oak_f_EC <- lmer(k ~ Location + (1|Replicate), data = Oak_f_EC)
summary(Lmer_oak_f_EC)
Anova(Lmer_oak_f_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_EC)
#qqplot
qqnorm(resid(Lmer_oak_f_EC))
qqline(resid(Lmer_oak_f_EC))
# transform k values
Oak_f_EC$Kt <- log10(Oak_f_EC$k)
# refit LMM
Lmer_oak_f_EC <- lmer(Kt ~ Location + (1|Replicate), data = Oak_f_EC, na.action = na.omit)
summary(Lmer_oak_f_EC)
Anova(Lmer_oak_f_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_EC)
#qqplot
qqnorm(resid(Lmer_oak_f_EC))
qqline(resid(Lmer_oak_f_EC))
#explained variation by model
r.squaredGLMM(Lmer_oak_f_EC)

#Oak Coarse EC
Oak_c_EC <- filter(Oak_EC, Mesh == "Coarse")
Lmer_oak_c_EC <- lmer(k ~ Location + (1|Replicate), data = Oak_c_EC)
summary(Lmer_oak_c_EC)
Anova(Lmer_oak_c_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_EC)
#qqplot
qqnorm(resid(Lmer_oak_c_EC))
qqline(resid(Lmer_oak_c_EC))
# transform k values
Oak_c_EC$Kt <- log10(Oak_c_EC$k)
# refit LMM
Lmer_oak_c_EC <- lmer(Kt ~ Location + (1|Replicate), data = Oak_c_EC, na.action = na.omit)
summary(Lmer_oak_c_EC)
Anova(Lmer_oak_c_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_EC)
#qqplot
qqnorm(resid(Lmer_oak_c_EC))
qqline(resid(Lmer_oak_c_EC))
#explained variation by model
r.squaredGLMM(Lmer_oak_c_EC)

#Wettex EC
Lmer_wettex_EC <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Wettex_EC)
summary(Lmer_wettex_EC)
Anova(Lmer_wettex_EC, test.statistic = "F")
#check model assumptions
plot(Lmer_wettex_EC)
#qqplot
qqnorm(resid(Lmer_wettex_EC))
qqline(resid(Lmer_wettex_EC))
# transform k values
Wettex_EC$Kt <- log10(Wettex_EC$k)
# refit LMM
Lmer_wettex_EC <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Wettex_EC, na.action = na.omit)
summary(Lmer_wettex_EC)
Anova(Lmer_wettex_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_EC)
#qqplot
qqnorm(resid(Lmer_wettex_EC))
qqline(resid(Lmer_wettex_EC))
#explained variation by model
r.squaredGLMM(Lmer_wettex_EC)

#Wettex fine EC
Wettex_f_EC <- filter(Wettex_EC, Mesh == "Fine")
Lmer_wettex_f_EC <- lmer(k ~ Location + (1|Replicate), data = Wettex_f_EC)
summary(Lmer_wettex_f_EC)
Anova(Lmer_wettex_f_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_EC)
#qqplot
qqnorm(resid(Lmer_wettex_f_EC))
qqline(resid(Lmer_wettex_f_EC))
# transform k values
Wettex_f_EC$Kt <- log10(Wettex_f_EC$k)
# refit LMM
Lmer_wettex_f_EC <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_f_EC, na.action = na.omit)
summary(Lmer_wettex_f_EC)
Anova(Lmer_wettex_f_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_EC)
#qqplot
qqnorm(resid(Lmer_wettex_f_EC))
qqline(resid(Lmer_wettex_f_EC))
#explained variation by model
r.squaredGLMM(Lmer_wettex_f_EC)

#Wettex coarse EC
Wettex_c_EC <- filter(Wettex_EC, Mesh == "Coarse")
Lmer_wettex_c_EC <- lmer(k ~ Location + (1|Replicate), data = Wettex_c_EC)
summary(Lmer_wettex_c_EC)
Anova(Lmer_wettex_c_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_EC)
#qqplot
qqnorm(resid(Lmer_wettex_c_EC))
qqline(resid(Lmer_wettex_c_EC))
# transform k values
Wettex_c_EC$Kt <- log10(Wettex_c_EC$k)
# refit LMM
Lmer_wettex_c_EC <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_c_EC, na.action = na.omit)
summary(Lmer_wettex_c_EC)
Anova(Lmer_wettex_c_EC, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_EC)
#qqplot
qqnorm(resid(Lmer_wettex_c_EC))
qqline(resid(Lmer_wettex_c_EC))
#explained variation by model
r.squaredGLMM(Lmer_wettex_c_EC)




##ENCLOSURE VS. DISTANCE
#Create separate dataframes for Oak and Wettex, Enclosure and Downstream
Oak_D <- filter(Bioassay_Data, Location == "Downstream", Bioassay == "Oak")
Wettex_D <- filter(Bioassay_Data, Location == "Downstream", Bioassay == "Wettex")
Oak_E <- filter(Bioassay_Data, Location == "Enclosure", Bioassay == "Oak")
Wettex_E <- filter(Bioassay_Data, Location == "Enclosure", Bioassay == "Wettex")

#Add distance
dist <- c(3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2,3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2)
Oak_D$Distance <- dist
Wettex_D$Distance <- dist
Oak_E$Distance <- 0
Wettex_E$Distance <- 0

#Combine Enclosure and Downstream
Oak_ED <- bind_rows(Oak_D, Oak_E)
Wettex_ED <- bind_rows(Wettex_D, Wettex_E)

#LMer Enclosure and Distance (EDis)
#Oak EDis
Lmer_oak_EDis <- lmer(k ~ Distance*Mesh + (1|Replicate), data = Oak_ED)
summary(Lmer_oak_EDis)
Anova(Lmer_oak_EDis, test.statistic = "F")
#check model assumptions
plot(Lmer_oak_ED)
#qqplot
qqnorm(resid(Lmer_oak_EDis))
qqline(resid(Lmer_oak_EDis))
# transform k values
Oak_ED$Kt <- log10(Oak_ED$k)
# refit LMM
Lmer_oak_EDis <- lmer(Kt ~ Distance*Mesh + (1|Replicate), data = Oak_ED, na.action = na.omit)
summary(Lmer_oak_EDis)
Anova(Lmer_oak_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_EDis)
#qqplot
qqnorm(resid(Lmer_oak_EDis))
qqline(resid(Lmer_oak_EDis))
# random-slope and random-intercept model
Lmer_oak_EDis_slope <- lmer(Kt ~ Distance*Mesh + (1 + Distance|Replicate), data = Oak_ED)
summary(Lmer_oak_EDis_slope)
anova(Lmer_oak_EDis, Lmer_oak_EDis_slope)
# Intercept-only model has lower AIC
#explained variation by model
r.squaredGLMM(Lmer_oak_EDis)

#Oak Fine EDis
Oak_f_EDis <- filter(Oak_ED, Mesh == "Fine")
Lmer_oak_f_EDis <- lmer(k ~ Distance + (1|Replicate), data = Oak_f_EDis)
summary(Lmer_oak_f_EDis)
Anova(Lmer_oak_f_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_EDis)
#qqplot
qqnorm(resid(Lmer_oak_f_EDis))
qqline(resid(Lmer_oak_f_EDis))
# transform k values
Oak_f_EDis$Kt <- log10(Oak_f_EDis$k)
# refit LMM
Lmer_oak_f_EDis <- lmer(Kt ~ Distance + (1|Replicate), data = Oak_f_EDis, na.action = na.omit)
summary(Lmer_oak_f_EDis)
Anova(Lmer_oak_f_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_EDis)
#qqplot
qqnorm(resid(Lmer_oak_f_EDis))
qqline(resid(Lmer_oak_f_EDis))
#explained variation by model
r.squaredGLMM(Lmer_oak_f_EDis)

#Oak Coarse EDis
Oak_c_EDis <- filter(Oak_ED, Mesh == "Coarse")
Lmer_oak_c_EDis <- lmer(k ~ Distance + (1|Replicate), data = Oak_c_EDis)
summary(Lmer_oak_c_EDis)
Anova(Lmer_oak_c_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_EDis)
#qqplot
qqnorm(resid(Lmer_oak_c_EDis))
qqline(resid(Lmer_oak_c_EDis))
# transform k values
Oak_c_EDis$Kt <- log10(Oak_c_EDis$k)
# refit LMM
Lmer_oak_c_EDis <- lmer(Kt ~ Distance + (1|Replicate), data = Oak_c_EDis, na.action = na.omit)
summary(Lmer_oak_c_EDis)
Anova(Lmer_oak_c_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_EDis)
#qqplot
qqnorm(resid(Lmer_oak_c_EDis))
qqline(resid(Lmer_oak_c_EDis))
#explained variation by model
r.squaredGLMM(Lmer_oak_c_EDis)

#Wettex EDis
Lmer_wettex_EDis <- lmer(k ~ Distance*Mesh + (1|Replicate), data = Wettex_ED)
summary(Lmer_wettex_EDis)
Anova(Lmer_wettex_EDis, test.statistic = "F")
#check model assumptions
plot(Lmer_wettex_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_EDis))
qqline(resid(Lmer_wettex_EDis))
# transform k values
Wettex_ED$Kt <- log10(Wettex_ED$k)
# refit LMM
Lmer_wettex_EDis <- lmer(Kt ~ Distance*Mesh + (1|Replicate), data = Wettex_ED, na.action = na.omit)
summary(Lmer_wettex_EDis)
Anova(Lmer_wettex_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_EDis))
qqline(resid(Lmer_wettex_EDis))
# random-slope and random-intercept model
Lmer_wettex_EDis_slope <- lmer(Kt ~ Distance*Mesh + (1 + Distance|Replicate), data = Wettex_ED)
summary(Lmer_wettex_EDis_slope)
anova(Lmer_wettex_EDis, Lmer_wettex_EDis_slope)
# Intercept-only model has lower AIC
#explained variation by model
r.squaredGLMM(Lmer_wettex_EDis)

#Wettex Fine EDis
Wettex_f_EDis <- filter(Wettex_ED, Mesh == "Fine")
Lmer_wettex_f_EDis <- lmer(k ~ Distance + (1|Replicate), data = Wettex_f_EDis)
summary(Lmer_wettex_f_EDis)
Anova(Lmer_wettex_f_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_f_EDis))
qqline(resid(Lmer_wettex_f_EDis))
# transform k values
Wettex_f_EDis$Kt <- log10(Wettex_f_EDis$k)
# refit LMM
Lmer_wettex_f_EDis <- lmer(Kt ~ Distance + (1|Replicate), data = Wettex_f_EDis, na.action = na.omit)
summary(Lmer_wettex_f_EDis)
Anova(Lmer_wettex_f_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_f_EDis))
qqline(resid(Lmer_wettex_f_EDis))
#explained variation by model
r.squaredGLMM(Lmer_wettex_f_EDis)

#Wettex Coarse EDis
Wettex_c_EDis <- filter(Wettex_ED, Mesh == "Coarse")
Lmer_wettex_c_EDis <- lmer(k ~ Distance + (1|Replicate), data = Wettex_c_EDis)
summary(Lmer_wettex_c_EDis)
Anova(Lmer_wettex_c_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_c_EDis))
qqline(resid(Lmer_wettex_c_EDis))
# transform k values
Wettex_c_EDis$Kt <- log10(Wettex_c_EDis$k)
# refit LMM
Lmer_wettex_c_EDis <- lmer(Kt ~ Distance + (1|Replicate), data = Wettex_c_EDis, na.action = na.omit)
summary(Lmer_wettex_c_EDis)
Anova(Lmer_wettex_c_EDis, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_EDis)
#qqplot
qqnorm(resid(Lmer_wettex_c_EDis))
qqline(resid(Lmer_wettex_c_EDis))
#explained variation by model
r.squaredGLMM(Lmer_wettex_c_EDis)





##ENCLOSURE VS. DOWNSTREAM
#LMer for Enclosure and Downstream (ED)
#Oak ED
Lmer_oak_ED <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Oak_ED)
summary(Lmer_oak_ED)
Anova(Lmer_oak_ED, test.statistic = "F")
#check model assumptions
plot(Lmer_oak_ED)
#qqplot
qqnorm(resid(Lmer_oak_ED))
qqline(resid(Lmer_oak_ED))
# transform k values
Oak_ED$Kt <- log10(Oak_ED$k)
# refit LMM
Lmer_oak_ED <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Oak_ED, na.action = na.omit)
summary(Lmer_oak_ED)
Anova(Lmer_oak_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_ED)
#qqplot
qqnorm(resid(Lmer_oak_ED))
qqline(resid(Lmer_oak_ED))
#explained variation by model
r.squaredGLMM(Lmer_oak_ED)

#Oak Fine ED 
Oak_f_ED <- filter(Oak_ED, Mesh == "Fine")
Lmer_oak_f_ED <- lmer(k ~ Location + (1|Replicate), data = Oak_f_ED)
summary(Lmer_oak_f_ED)
Anova(Lmer_oak_f_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_ED)
#qqplot
qqnorm(resid(Lmer_oak_f_ED))
qqline(resid(Lmer_oak_f_ED))
# transform k values
Oak_f_ED$Kt <- log10(Oak_f_ED$k)
# refit LMM
Lmer_oak_f_ED <- lmer(Kt ~ Location + (1|Replicate), data = Oak_f_ED, na.action = na.omit)
summary(Lmer_oak_f_ED)
Anova(Lmer_oak_f_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_ED)
#qqplot
qqnorm(resid(Lmer_oak_f_ED))
qqline(resid(Lmer_oak_f_ED))
#explained variation by model
r.squaredGLMM(Lmer_oak_f_ED)

#Oak Coarse ED
Oak_c_ED <- filter(Oak_ED, Mesh == "Coarse")
Lmer_oak_c_ED <- lmer(k ~ Location + (1|Replicate), data = Oak_c_ED)
summary(Lmer_oak_c_ED)
Anova(Lmer_oak_c_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_ED)
#qqplot
qqnorm(resid(Lmer_oak_c_ED))
qqline(resid(Lmer_oak_c_ED))
# transform k values
Oak_c_ED$Kt <- log10(Oak_c_ED$k)
# refit LMM
Lmer_oak_c_ED <- lmer(Kt ~ Location + (1|Replicate), data = Oak_c_ED, na.action = na.omit)
summary(Lmer_oak_c_ED)
Anova(Lmer_oak_c_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_ED)
#qqplot
qqnorm(resid(Lmer_oak_c_ED))
qqline(resid(Lmer_oak_c_ED))
#explained variation by model
r.squaredGLMM(Lmer_oak_c_ED)

#Wettex ED
Lmer_wettex_ED <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Wettex_ED)
summary(Lmer_wettex_ED)
Anova(Lmer_wettex_ED, test.statistic = "F")
#check model assumptions
plot(Lmer_wettex_ED)
#qqplot
qqnorm(resid(Lmer_wettex_ED))
qqline(resid(Lmer_wettex_ED))
# transform k values
Wettex_ED$Kt <- log10(Wettex_ED$k)
# refit LMM
Lmer_wettex_ED <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Wettex_ED, na.action = na.omit)
summary(Lmer_wettex_ED)
Anova(Lmer_wettex_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_ED)
#qqplot
qqnorm(resid(Lmer_wettex_ED))
qqline(resid(Lmer_wettex_ED))
#explained variation by model
r.squaredGLMM(Lmer_wettex_ED)

#Wettex fine ED
Wettex_f_ED <- filter(Wettex_ED, Mesh == "Fine")
Lmer_wettex_f_ED <- lmer(k ~ Location + (1|Replicate), data = Wettex_f_ED)
summary(Lmer_wettex_f_ED)
Anova(Lmer_wettex_f_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_ED)
#qqplot
qqnorm(resid(Lmer_wettex_f_ED))
qqline(resid(Lmer_wettex_f_ED))
# transform k values
Wettex_f_ED$Kt <- log10(Wettex_f_ED$k)
# refit LMM
Lmer_wettex_f_ED <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_f_ED, na.action = na.omit)
summary(Lmer_wettex_f_ED)
Anova(Lmer_wettex_f_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_ED)
#qqplot
qqnorm(resid(Lmer_wettex_f_ED))
qqline(resid(Lmer_wettex_f_ED))
#explained variation by model
r.squaredGLMM(Lmer_wettex_f_ED)

#Wettex coarse ED
Wettex_c_ED <- filter(Wettex_ED, Mesh == "Coarse")
Lmer_wettex_c_ED <- lmer(k ~ Location + (1|Replicate), data = Wettex_c_ED)
summary(Lmer_wettex_c_ED)
Anova(Lmer_wettex_c_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_ED)
#qqplot
qqnorm(resid(Lmer_wettex_c_ED))
qqline(resid(Lmer_wettex_c_ED))
# transform k values
Wettex_c_ED$Kt <- log10(Wettex_c_ED$k)
# refit LMM
Lmer_wettex_c_ED <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_c_ED, na.action = na.omit)
summary(Lmer_wettex_c_ED)
Anova(Lmer_wettex_c_ED, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_ED)
#qqplot
qqnorm(resid(Lmer_wettex_c_ED))
qqline(resid(Lmer_wettex_c_ED))
#explained variation by model
r.squaredGLMM(Lmer_wettex_c_ED)





##CONTROL VS. DOWNSTREAM
#Create separate dataframes for Oak and Wettex, Control and Downstream
Bioassay_Data_CD <- dplyr::filter(Bioassay_Data, Location == "Control" | Location == "Downstream")
Oak_CD <- dplyr::filter(Bioassay_Data_CD, Bioassay == "Oak")
Wettex_CD <- dplyr::filter(Bioassay_Data_CD, Bioassay == "Wettex") 

#LMM control vs. downstream (CD)
#Oak CD
Lmer_oak_CD <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Oak_CD)
summary(Lmer_oak_CD)
Anova(Lmer_oak_CD, test.statistic = "F")
#check model assumptions
plot(Lmer_oak_CD)
#qqplot
qqnorm(resid(Lmer_oak_CD))
qqline(resid(Lmer_oak_CD))
# transform k values
Oak_CD$Kt <- log10(Oak_CD$k)
# refit LMM
Lmer_oak_CD <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Oak_CD, na.action = na.omit)
summary(Lmer_oak_CD)
Anova(Lmer_oak_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_CD)
#qqplot
qqnorm(resid(Lmer_oak_CD))
qqline(resid(Lmer_oak_CD))
#explained variation by model
r.squaredGLMM(Lmer_oak_CD)

#Oak Fine CD
Oak_f_CD <- filter(Oak_CD, Mesh == "Fine")
Lmer_oak_f_CD <- lmer(k ~ Location + (1|Replicate), data = Oak_f_CD)
summary(Lmer_oak_f_CD)
Anova(Lmer_oak_f_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_CD)
#qqplot
qqnorm(resid(Lmer_oak_f_CD))
qqline(resid(Lmer_oak_f_CD))
# transform k values
Oak_f_CD$Kt <- log10(Oak_f_CD$k)
# refit LMM
Lmer_oak_f_CD <- lmer(Kt ~ Location + (1|Replicate), data = Oak_f_CD, na.action = na.omit)
summary(Lmer_oak_f_CD)
Anova(Lmer_oak_f_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_f_CD)
#qqplot
qqnorm(resid(Lmer_oak_f_CD))
qqline(resid(Lmer_oak_f_CD))
#explained variation by model
r.squaredGLMM(Lmer_oak_f_CD)

#Oak Coarse CD
Oak_c_CD <- filter(Oak_CD, Mesh == "Coarse")
Lmer_oak_c_CD <- lmer(k ~ Location + (1|Replicate), data = Oak_c_CD)
summary(Lmer_oak_c_CD)
Anova(Lmer_oak_c_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_CD)
#qqplot
qqnorm(resid(Lmer_oak_c_CD))
qqline(resid(Lmer_oak_c_CD))
# transform k values
Oak_c_CD$Kt <- log10(Oak_c_CD$k)
# refit LMM
Lmer_oak_c_CD <- lmer(Kt ~ Location + (1|Replicate), data = Oak_c_CD, na.action = na.omit)
summary(Lmer_oak_c_CD)
Anova(Lmer_oak_c_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_oak_c_CD)
#qqplot
qqnorm(resid(Lmer_oak_c_CD))
qqline(resid(Lmer_oak_c_CD))
#explained variation by model
r.squaredGLMM(Lmer_oak_c_CD)

#Wettex CD
Lmer_wettex_CD <- lmer(k ~ Location*Mesh + (1|Site/Replicate), data = Wettex_CD)
summary(Lmer_wettex_CD)
Anova(Lmer_wettex_CD, test.statistic = "F")
#check model assumptions
plot(Lmer_wettex_CD)
#qqplot
qqnorm(resid(Lmer_wettex_CD))
qqline(resid(Lmer_wettex_CD))
# transform k values
Wettex_CD$Kt <- log10(Wettex_CD$k)
# refit LMM
Lmer_wettex_CD <- lmer(Kt ~ Location*Mesh + (1|Site/Replicate), data = Wettex_CD, na.action = na.omit)
summary(Lmer_wettex_CD)
Anova(Lmer_wettex_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_CD)
#qqplot
qqnorm(resid(Lmer_wettex_CD))
qqline(resid(Lmer_wettex_CD))
#explained variation by model
r.squaredGLMM(Lmer_wettex_CD)

#Wettex fine CD
Wettex_f_CD <- filter(Wettex_CD, Mesh == "Fine")
Lmer_wettex_f_CD <- lmer(k ~ Location + (1|Replicate), data = Wettex_f_CD)
summary(Lmer_wettex_f_CD)
Anova(Lmer_wettex_f_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_CD)
#qqplot
qqnorm(resid(Lmer_wettex_f_CD))
qqline(resid(Lmer_wettex_f_CD))
# transform k values
Wettex_f_CD$Kt <- log10(Wettex_f_CD$k)
# refit LMM
Lmer_wettex_f_CD <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_f_CD, na.action = na.omit)
summary(Lmer_wettex_f_CD)
Anova(Lmer_wettex_f_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_f_CD)
#qqplot
qqnorm(resid(Lmer_wettex_f_CD))
qqline(resid(Lmer_wettex_f_CD))
#explained variation by model
r.squaredGLMM(Lmer_wettex_f_CD)

#Wettex coarse CD
Wettex_c_CD <- filter(Wettex_CD, Mesh == "Coarse")
Lmer_wettex_c_CD <- lmer(k ~ Location + (1|Replicate), data = Wettex_c_CD)
summary(Lmer_wettex_c_CD)
Anova(Lmer_wettex_c_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_CD)
#qqplot
qqnorm(resid(Lmer_wettex_c_CD))
qqline(resid(Lmer_wettex_c_CD))
# transform k values
Wettex_c_CD$Kt <- log10(Wettex_c_CD$k)
# refit LMM
Lmer_wettex_c_CD <- lmer(Kt ~ Location + (1|Replicate), data = Wettex_c_CD, na.action = na.omit)
summary(Lmer_wettex_c_CD)
Anova(Lmer_wettex_c_CD, test.statistic = "F")
# check model assumptions
plot(Lmer_wettex_c_CD)
#qqplot
qqnorm(resid(Lmer_wettex_c_CD))
qqline(resid(Lmer_wettex_c_CD))
#explained variation by model
r.squaredGLMM(Lmer_wettex_c_CD)






##GREEN TEA ANALYSIS
#Load Data
Green_Tea_Data <- read.csv(file = "Green Tea Data.csv")

#Format Data
tea_k<-as.data.frame(Green_Tea_Data)
tea_k$Site<-as.factor(tea_k$Site)

#Simplify df and pivot
tea_k <- dplyr::select(tea_k, Site, Location, type, k_R1, k_R2, k_R3, k_R4) %>% 
  dplyr::rename(R1 = k_R1,
                R2 = k_R2,
                R3 = k_R3,
                R4 = k_R4,
  ) %>% 
  pivot_longer(
    cols = starts_with("R"), 
    names_to = "Replicate", 
    values_to = "k",
    values_drop_na = TRUE
  )

#ENCLOSURE VS. CONTROL
#Separate data frame for Enclosure and Control
tea_EC <- dplyr::filter(tea_k, Location == "Enclosure" | Location == "Control")

#Enclosure vs. Control LMM
hist(tea_EC$k)
lmer_tea_EC <- lmer(k ~ Location + (1|Site) + (1|Replicate), data = tea_EC)
summary(lmer_tea_EC)
Anova(lmer_tea_EC, test.statistic = "F")
# check model assumptions
plot(lmer_tea_EC)
#qqplot
qqnorm(resid(lmer_tea_EC))
qqline(resid(lmer_tea_EC))
#explained variation by model
r.squaredGLMM(lmer_tea_EC)
# random-slope and random-intercept model
lmer_tea_EC_Slope <- lmer(k ~ Location + (1 + Location|Site) + (1 + Location|Replicate), data = tea_EC)
summary(lmer_tea_EC_Slope)
anova(lmer_tea_EC, lmer_tea_EC_Slope)
# no significant difference between models
AIC(lmer_tea_EC, lmer_tea_EC_Slope, k = T)
# intercept-only model has lower AIC
# pairwise
lmer_tea_sites <- lmer(k ~ Site + (1|Replicate), data = tea_EC)
em <- emmeans(lmer_tea_sites, c("Site"))
contrast(em, method = "pairwise")


#ENCLOSURE VS. DISTANCE
#Separate data frames for enclosure & downstream
tea_D <- filter(tea_k, Location == "Downstream")
tea_E <- filter(tea_k, Location == "Enclosure")

#add distance
dist <- c(3.2, 3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.2)
tea_D$Distance <-dist
tea_E$Distance <- 0

#combine D and E
tea_ED <- bind_rows(tea_D, tea_E)

#Enclosure and Distance LMM
lm_EDis<-lmer(k ~ Distance + (1|Replicate), data = tea_ED)
summary(lm_EDis)
Anova(lm_EDis, test.statistic = "F")
# check model assumptions
plot(lm_EDis)
summary(lm_EDis)
#explained variation by model
r.squaredGLMM(lm_EDis)
# random-slope and random-intercept model
lm_EDis_Slope <- lmer(k ~ Distance + (1 + Distance|Replicate), data = tea_ED)
summary(lm_EDis_Slope)
anova(lm_EDis, lm_EDis_Slope)
# no significant difference between models
AIC(lm_EDis, lm_EDis_Slope, k = 2)
# AIC are equal


#ENCLOSURE VS. DOWNSTREAM
#Enclosure and Downstream LMM
hist(tea_ED$k)
lmer_tea_ED <- lmer(k ~ Location + (1|Site) + (1|Replicate), data = tea_ED)
summary(lmer_tea_ED)
Anova(lmer_tea_ED, test.statistic = "F")
# check model assumptions
plot(lmer_tea_ED)
#qqplot
qqnorm(resid(lmer_tea_ED))
qqline(resid(lmer_tea_ED))
#explained variation by model
r.squaredGLMM(lmer_tea_ED)

#Separate data frame for Control and Downstream
tea_CD <- dplyr::filter(tea_k, Location == "Control" | Location == "Downstream")

#Control and Downstream LMM
hist(tea_CD$k)
lmer_tea_CD <- lmer(k ~ Location + (1|Site) + (1|Replicate), data = tea_CD)
summary(lmer_tea_CD)
Anova(lmer_tea_CD, test.statistic = "F")
# check model assumptions
plot(lmer_tea_CD)
#qqplot
qqnorm(resid(lmer_tea_CD))
qqline(resid(lmer_tea_CD))
#explained variation by model
r.squaredGLMM(lmer_tea_CD)





##SCRIPT FOR GRAPHS/FIGURES
#Barcharts for Oak and Wettex/Enclosure vs. Control
head(Bioassay_EC_mean_SE)
Bioassay_EC_bar <- ggplot(Bioassay_EC_mean_SE, aes(x=Location, y=mean_k, color=Location)) +
  geom_bar(stat="identity", fill="skyblue", alpha = 0.7) +
  geom_errorbar( aes(x=Location, ymin=mean_k-SE_k, ymax=mean_k+SE_k), width = 0.4, colour="black", alpha=0.9, size=1.3) +
  facet_grid(. ~ Mesh + Bioassay) +
  labs(y = expression(k[combined]), title = "\n Enclosure vs. Control \n") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
Bioassay_EC_bar

#Barcharts for Oak and Wettex/Enclosure vs. Downstream
head(Bioassay_ED_mean_SE)
Bioassay_ED_bar <- ggplot(Bioassay_ED_mean_SE, aes(x=Location, y=mean_k, color=Location)) +
  geom_bar(stat="identity", fill="skyblue", alpha = 0.7) +
  geom_errorbar( aes(x=Location, ymin=mean_k-SE_k, ymax=mean_k+SE_k), width = 0.4, colour="black", alpha=0.9, size=1.3) +
  facet_grid(. ~ Mesh + Bioassay) +
  labs(y = expression(k[combined]), title = "\n Enclosure vs. Downstream \n") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
Bioassay_ED_bar

#Barcharts for Oak and Wettex/Control vs. Downstream
head(Bioassay_CD_mean_SE)
Bioassay_CD_bar <- ggplot(Bioassay_CD_mean_SE, aes(x=Location, y=mean_k, color=Location)) +
  geom_bar(stat="identity", fill="skyblue", alpha = 0.7) +
  geom_errorbar( aes(x=Location, ymin=mean_k-SE_k, ymax=mean_k+SE_k), width = 0.4, colour="black", alpha=0.9, size=1.3) +
  facet_grid(. ~ Mesh + Bioassay) +
  labs(y = expression(k[combined]), title = "\n Control vs. Downstream \n") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
Bioassay_CD_bar

#Barcharts for microbial breakdown only (including green tea)
head(Bioassay_microbial)
Bioassay_microbial_bar <- ggplot(Bioassay_microbial, aes(x=Location, y=mean_k, color=Location)) +
  geom_bar(stat="identity", fill="skyblue", alpha = 0.7) +
  geom_errorbar( aes(x=Location, ymin=mean_k-SE_k, ymax=mean_k+SE_k), width = 0.4, colour="black", alpha=0.9, size=1.3) +
  facet_grid(. ~ C_N + Bioassay) +
  labs(x = "Treatment", y = "k-microbial", title = "\n Microbial breakdown \n") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
Bioassay_microbial_bar

#Scatterplot for microbial breakdown 
C_N_data$C_N <- as.factor(C_N_data$C_N)
C_N_data$Bioassay <- as.factor(C_N_data$Bioassay)
C_N_line <- ggplot(C_N_data, aes(x=C_N, y=k)) +
  geom_point(aes(color=Bioassay)) +
  geom_smooth(method = "lm",
              se = TRUE,
              color="black") +
  labs(title = "Bioassay C:N ratio against Microbial breakdown rate",
       x = "C:N ratio",
       y = "Microbial breakdown rate")
C_N_line

#Combine Enclosure and Downstream data for Oak and Wettex
Bioassay_Data_EDis <- bind_rows(Oak_ED, Wettex_ED)
#Filter for microbial breakdown only
Bioassay_f_EDis <- filter(Bioassay_Data_EDis, Mesh == "Fine")
#rename columns in green tea data
dat <- data.frame(old=c("type"), new=c("Bioassay"), stringsAsFactors=FALSE)
match(dat[,"old"], names(tea_ED))
names(tea_ED)[match(dat[,"old"], names(tea_ED))] = dat[,"new"]

#Combine Green tea data to Oak and Wettex
Bioassay_f_EDis <- bind_rows(Bioassay_f_EDis, tea_ED)
Bioassay_f_Edis <- filter(Bioassay_f_EDis, Bioassay == "Oak" | Bioassay == "Green")

#Scatterplot for Enclosure vs. Distance
Edis_line <- ggplot(Bioassay_f_Edis, aes(x=Distance, y=k)) +
  geom_point(aes(shape = Bioassay)) +
  geom_point(aes(color = Bioassay)) +
  geom_smooth(method = "lm",
              se = TRUE,
              aes(color = factor(Bioassay))) +
  labs(title = "Enclosure vs Distance",
       x = "Distance (km)",
       y = expression(k[microbial]))
Edis_line <- Edis_line + scale_fill_discrete(name="Bioassay")
Edis_line



##INVERTEBRATE ANALYSIS (GAMMARUS ABUNDANCE)
#Load Data
Inverts <- read.csv(file = "Gammarus abundance.csv")

#Format Data
Inverts$Site <- as.factor(Inverts$Site)
Inverts$Bioassay <- as.factor(Inverts$Bioassay)
Inverts <- dplyr::select(Inverts, Site, Location, Bioassay, Replicate, k, Adult_gammarus, Juvenile_gammarus, Xt)

# new column with correction to give abundance per unit mass of leaf litter remaining in bag
Inverts$abundance_adult_gammarus <- Inverts$Adult_gammarus/Inverts$Xt
Inverts$abundance_juvenile_gammarus <- Inverts$Juvenile_gammarus/Inverts$Xt

# remove 0s
Inverts[Inverts==0] <- NA

#Enclosure and Control
Inverts_Oak <- filter(Inverts, Bioassay == "Oak", Location == "Enclosure" | Location == "Control")
Inverts_Wettex <- filter(Inverts, Bioassay == "Wettex", Location == "Enclosure" | Location == "Control")
Inverts_Bioassay <- filter(Inverts, Location == "Enclosure" | Location == "Control")

#Adult gammarus
lmer_adult_gammarus <- lmer(k~abundance_adult_gammarus*Bioassay*Location + (1|Site/Replicate), data = Inverts_Bioassay, na.action = na.omit)
summary(lmer_adult_gammarus)
Anova(lmer_adult_gammarus)
# check model assumptions
plot(lmer_adult_gammarus)
#qqplot
qqnorm(resid(lmer_adult_gammarus))
qqline(resid(lmer_adult_gammarus))
#explained variation by model
r.squaredGLMM(lmer_adult_gammarus)

#Juvenile gammarus
lmer_juvenile_gammarus <- lmer(k~abundance_juvenile_gammarus*Bioassay*Location + (1|Site/Replicate), data = Inverts_Bioassay, na.action = na.omit)
summary(lmer_juvenile_gammarus)
Anova(lmer_juvenile_gammarus)
# check model assumptions
plot(lmer_juvenile_gammarus)
#qqplot
qqnorm(resid(lmer_juvenile_gammarus))
qqline(resid(lmer_juvenile_gammarus))
#explained variation by model
r.squaredGLMM(lmer_juvenile_gammarus)

#Adult gammarus Oak
lmer_adult_gammarus_oak <- lmer(k~abundance_adult_gammarus*Location + (1|Replicate), data = Inverts_Oak, na.action = na.omit)
summary(lmer_adult_gammarus_oak)
Anova(lmer_adult_gammarus_oak)
# check model assumptions
plot(lmer_adult_gammarus_oak)
#qqplot
qqnorm(resid(lmer_adult_gammarus_oak))
qqline(resid(lmer_adult_gammarus_oak))
#explained variation by model
r.squaredGLMM(lmer_adult_gammarus_oak)

#Adult gammarus Wettex
lmer_adult_gammarus_wettex <- lmer(k~abundance_adult_gammarus*Location + (1|Replicate), data = Inverts_Wettex, na.action = na.omit)
summary(lmer_adult_gammarus_wettex)
Anova(lmer_adult_gammarus_wettex)
# check model assumptions
plot(lmer_adult_gammarus_wettex)
#qqplot
qqnorm(resid(lmer_adult_gammarus_wettex))
qqline(resid(lmer_adult_gammarus_wettex))
#explained variation by model
r.squaredGLMM(lmer_adult_gammarus_wettex)

#Juvenile gammarus Oak
lmer_juvenile_gammarus_oak <- lmer(k~abundance_juvenile_gammarus*Location + (1|Replicate), data = Inverts_Oak, na.action = na.omit)
summary(lmer_juvenile_gammarus_oak)
Anova(lmer_juvenile_gammarus_oak)
# check model assumptions
plot(lmer_juvenile_gammarus_oak)
#qqplot
qqnorm(resid(lmer_juvenile_gammarus_oak))
qqline(resid(lmer_juvenile_gammarus_oak))
#explained variation by model
r.squaredGLMM(lmer_juvenile_gammarus_oak)

#Juvenile gammarus Wettex
lmer_juvenile_gammarus_wettex <- lmer(k~abundance_juvenile_gammarus*Location + (1|Replicate), data = Inverts_Wettex, na.action = na.omit)
summary(lmer_juvenile_gammarus_wettex)
Anova(lmer_juvenile_gammarus_wettex)
# check model assumptions
plot(lmer_juvenile_gammarus_wettex)
#qqplot
qqnorm(resid(lmer_juvenile_gammarus_wettex))
qqline(resid(lmer_juvenile_gammarus_wettex))
#explained variation by model
r.squaredGLMM(lmer_juvenile_gammarus_wettex)

#Adult gammarus scatterplot with line of best fit
AGammarus_EC_line <- ggplot(Inverts_Bioassay, aes(x=abundance_adult_gammarus, y=k)) +
  geom_point(aes(shape = Bioassay)) +
  geom_point(aes(color = Bioassay)) +
  geom_smooth(method = "lm",
              se = TRUE,
              aes(color = factor(Bioassay))) +
  labs(title = "Gammarus abundance against breakdown rate",
       x = "Adult gammarus abundance",
       y = expression(k[total]))
AGammarus_EC_line

#Juvenile gammarus scatterplot with line of best fit
JGammarus_EC_line <- ggplot(Inverts_Bioassay, aes(x=abundance_juvenile_gammarus, y=k)) +
  geom_point(aes(shape = factor(Bioassay))) +
  geom_point(aes(color = factor(Bioassay))) +
  geom_smooth(method = "lm",
              se = TRUE,
              aes(color = factor(Bioassay))) +
  labs(title = "Gammarus abundance against breakdown rate",
       x = "Juvenile gammarus abundance",
       y = "Total breakdown rate")
JGammarus_EC_line





##INVERTEBRATE PCA
#Load Data
invert_site_values <- read.csv(file = "invert pca.csv")

#Format Data
invert_df<-as.data.frame(invert_site_values)
invert_pca <- princomp(invert_df[,-c(1,10)], cor = TRUE)

#Evaluation of the Axes
#1. The Latent Root Criterion
eigenvalues <- invert_pca$sdev^2
eigenvalues
#Component 1,2 and 3 have eigenvalues above 1
#This Criterion suggest three components are sufficient

#2. The Scree Plot Criterion
plot(invert_pca, type="lines", ylim=c(0,3))
#The 'elbow' of the scree plot is at comp.3
#This criterion suggests three components are sufficient

#3. The Relative Percent Variance Criterion
summary(invert_pca)
#The first three components are sufficient in explaining 70.26% of variation
#All three criterion suggest three components are sufficient

#Interpreting of the Axes
loadings(invert_pca)

#gammarus - comp.1
#acellus - comp.2
#cased caddis - comp.2
#baby gammarus - comp.1
#caseless cadis - comp.3
#baby crayfish - comp.3
#leech - comp.3
#snails - comp.2

#You could group 
#Comp.1 as Primary shredders
#Comp.2 as Secondary shredders
#Comp.3 as Predators

#plot pca
biplot(invert_pca)
fviz_pca_biplot(invert_pca, repel = TRUE, geom.ind = "point", ellipse.level=0.95, col.var = "black", labelsize=4)

invert.treatment <- c(rep("Downstream", 5), rep("Enclosure", 5), rep("Control", 5))
invert_df <- cbind(invert_df, invert.treatment)
ggbiplot(invert_pca)
ggbiplot(invert_pca, ellipse=TRUE, groups=invert.treatment)
ggbiplot(invert_pca, ellipse=TRUE, choices=c(1,3), groups=invert.treatment)
ggbiplot(invert_pca, ellipse=TRUE, choices=c(2,3), groups=invert.treatment)

#pca summary stats
prin_comp_inv <- prcomp(invert_df[,-c(1,10)], rank. = 3)
components <- prin_comp_inv[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, invert_df$invert.treatment)
tot_explained_variance_ratio <- summary(prin_comp_inv)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

#attempt at 3D plot
tit = 'Total Explained Variance = 99.998'
Three_D_plot <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~invert_df$invert.treatment, colors = c('#636EFA','#EF553B','#00CC96') ) %>%
  add_markers(size = 12)
Three_D_plot <- Three_D_plot %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#e5ecf6")
  )
Three_D_plot




##BOD ANALYSIS

#BOD Anova
BOD_t <- as.data.frame(BOD_t)
BOD_t$Site <- as.factor(BOD_t$Site)
BOD_t$Treatment <- as.factor(BOD_t$Treatment)

#Simplify df and pivot
BOD_t<- dplyr::select(BOD_t, Site, Treatment, BOD5_t1, BOD5_t2, BOD5_t3) %>%
  dplyr::rename(R1 = BOD5_t1,
                R2 = BOD5_t2,
                R3 = BOD5_t3
  ) %>%
  pivot_longer(
    cols=starts_with("R"),
    names_to = "Replicate",
    values_to = "BOD5",
    values_drop_na = TRUE)

#Separate dataframes for analysis
BOD_t_EC <- dplyr::filter(BOD_t, Treatment == "Enclosure" | Treatment == "Control")
BOD_t_ED <- dplyr::filter(BOD_t, Treatment == "Enclosure" | Treatment == "Downstream")
BOD_t_CD <- dplyr::filter(BOD_t, Treatment == "Control" | Treatment == "Downstream")

#Create separate dataframes for Enclosure and Downstream
BOD_D <- filter(BOD_t, Treatment == "Downstream")
BOD_E <- filter(BOD_t, Treatment == "Enclosure")

#Add distance
Dist <- c(3.2, 3.2, 3.2, 1.6, 1.6, 1.6, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2)
BOD_D$Distance <- Dist
BOD_E$Distance <- 0

#Combine Enclosure and Downstream
BOD_ED <- bind_rows(BOD_D, BOD_E)

#BOD Enclosure vs. Control
Lmer_BOD_EC <- lmer(BOD5 ~ Treatment + (1|Replicate), data = BOD_t_EC)
summary(Lmer_BOD_EC)
Anova(Lmer_BOD_EC, test.statistic = "F")
#check model assumptions
plot(Lmer_BOD_EC)
#qqplot
qqnorm(resid(Lmer_BOD_EC))
qqline(resid(Lmer_BOD_EC))
#explained variation by model
r.squaredGLMM(Lmer_BOD_EC)

#BOD Enclosure vs. Downstream
Lmer_BOD_ED <- lmer(BOD5 ~ Treatment + (1|Replicate), data = BOD_t_ED)
summary(Lmer_BOD_ED)
Anova(Lmer_BOD_ED, test.statistic = "F")
#check model assumptions
plot(Lmer_BOD_ED)
#qqplot
qqnorm(resid(Lmer_BOD_ED))
qqline(resid(Lmer_BOD_ED))
#explained variation by model
r.squaredGLMM(Lmer_BOD_ED)

#BOD Enclosure vs. Distance
Lmer_BOD_EDis <- lmer(BOD5 ~ Distance + (1|Replicate), data = BOD_ED)
summary(Lmer_BOD_EDis)
Anova(Lmer_BOD_EDis, test.statistic = "F")
#check model assumptions
plot(Lmer_BOD_EDis)
#qqplot
qqnorm(resid(Lmer_BOD_EDis))
qqline(resid(Lmer_BOD_EDis))
#explained variation by model
r.squaredGLMM(Lmer_BOD_EDis)

#BOD Control vs. Downstream
Lmer_BOD_CD <- lmer(BOD5 ~ Treatment + (1|Replicate), data = BOD_t_CD)
summary(Lmer_BOD_CD)
Anova(Lmer_BOD_CD, test.statistic = "F")
#check model assumptions
plot(Lmer_BOD_CD)
#qqplot
qqnorm(resid(Lmer_BOD_CD))
qqline(resid(Lmer_BOD_CD))
#explained variation by model
r.squaredGLMM(Lmer_BOD_CD)

#Boxplot
BOD_df<-as.data.frame(BOD_t)
level_BOD <- c('Downstream', 'Enclosure', 'Control')
BOD_df$Treatment <- c('Downstream','Downstream','Downstream','Downstream','Downstream', 
                      'Enclosure','Enclosure','Enclosure','Enclosure','Enclosure', 'Control', 'Control', 'Control', 'Control', 'Control')

BOD_df_mean<-data.frame(Mean=c(3.19467,3.54333,3.13267),
                        se=c(0.0979,0.097,0.10542),
                        Treatment=as.factor(c("Downstream", "Enclosure", "Control")),
                        Category=c("Sites 1-5", "Sites 6-10", "Sites 11-15"))


ggplot(BOD_df_mean, aes(x=Category, y=Mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.2) +
  labs(title = "Biochemical oxygen demand",
       x = element_blank(),
       y = "mean BOD5")