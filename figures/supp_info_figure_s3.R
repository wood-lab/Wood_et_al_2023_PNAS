

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

data<-readRDS("data/compiled_data.RDS")


# Isolate the pollutants from the larger dataset

pollutants <- cbind(data$year,data$pb,data$as,data$zn,data$ni,data$v,data$cr,data$cu,data$ba,data$be,
                    data$sig8_lignin,data$lamb8,data$bd.v_soil_biomarker)
colnames(pollutants)<-c("year","Pb","As","Zn","Ni","V","Cr","Cu","Ba","Be","Sig8_lignin","Lamb8","Bd.V_soil_biomarker")


# Make correlation plot

matrix<-cor(pollutants,use="pairwise.complete.obs")
par(xpd = TRUE)
corrplot(matrix,method="number",type="lower",tl.col="black",tl.srt=45,mar=c(0,0,0,0))

