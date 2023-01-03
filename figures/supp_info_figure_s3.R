### Supplementary Information Figure S3
### Created by Chelsea Wood
### chelwood@uw.edu

library(tidyverse)
library(dplyr)
library(DHARMa)
library(lmtest)
library(glmmTMB)
library(reshape2)
library(gamm4)
library(MASS)
library(sjPlot)
library(car)
library(corrplot)
select <- dplyr::select

# Load pollutant data

core_dates<-read.csv("data/contaminants_brandenberger/core_dates.csv",header=T,sep=",")
core_dates$unique_ID<-paste(core_dates$Core_ID,core_dates$Core_section_depth,sep="_")
core_metals<-read.csv("data/contaminants_brandenberger/core_metals.csv",header=T,sep=",")
core_metals$unique_ID<-paste(core_metals$Core_ID,core_metals$Core_section,sep="_")
core_organics<-read.csv("data/contaminants_brandenberger//core_organics.csv",header=T,sep=",")
core_organics$unique_ID<-paste(core_organics$Core_ID,core_organics$Core_section,sep="_")
core_data_metals<-merge(core_dates,core_metals,by="unique_ID")
core_data_organics<-merge(core_dates,core_organics,by="unique_ID")

# Average the two cores

core_data_averaged_metals <- core_data_metals %>%
  group_by(Est_year) %>%
  summarize_at(vars(Pb,As,Zn,Ni,V,Cr,Cu,Ba,Be),funs(mean = mean(.,na.rm=TRUE)))

core_data_averaged_organics <- core_data_organics %>%
  group_by(Est_year) %>%
  summarize_at(vars(Sig8_lignin,Lamb8,Bd.V_soil_biomarker),funs(mean = mean(.,na.rm=TRUE)))

core_data_averaged <- merge(core_data_averaged_metals,core_data_averaged_organics,by="Est_year",all=TRUE)
colnames(core_data_averaged)<-c("year","Pb","As","Zn","Ni","V","Cr","Cu","Ba","Be","Sig8_lignin","Lamb8","Bd.V_soil_biomarker")


# Make correlation plot

matrix<-cor(core_data_averaged,use="pairwise.complete.obs")
par(xpd = TRUE)
corrplot(matrix,method="number",type="lower",tl.col="black",tl.srt=45,mar=c(0,0,0,0))

