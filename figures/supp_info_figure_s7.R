### Supplementary Information Figure S7
### Created by Tim Essington
### essing@uw.edu

# Script to make QQ plot for Model 3 (phase 1) and temp model (phase 2)
library(TMB)
library(glmmTMB)
library(tidyr)
library(ggplot2)
library(dplyr)
library(DHARMa)
source("src/DHARMa_funs.R")
# set simulation values
nSims <- 250

# load data
source("src/load_data.R")


### Diagnostics for phase 1 models ####
source("src/configure_TMB_runs.R")


##### Load model output #####
load("src/analysis/phase_2_results.Rdata")
model <-models$temp$tempyear
dat.list <- model$dat
rep <- model$rep

datnames <- names(dat.list)
for (i in 1:length(datnames)) assign(datnames[i], dat.list[[i]])

#### Extract Parameters ####

fixed_effect <- summary(model$rep, "fixed")
fixed_names <- rownames(fixed_effect)
rep_effect <- summary(model$rep, "report")
rep_names <- rownames(rep_effect)
ra_effect <- summary(model$rep, "random")
ra_names <- rownames(ra_effect)

# 
phi <- rep_effect[grepl(rep_names, pattern = "phi"),1]
bo <-  ra_effect[grepl(ra_names, pattern = "\\bbo\\b"),1]
blat<-  ra_effect[grepl(ra_names, pattern = "\\bblat\\b"),1]
blength<-  ra_effect[grepl(ra_names, pattern = "\\bblength\\b"),1]
gbar <- fixed_effect[grepl(fixed_names, pattern = "\\bgbar\\b"),1]
z <- ra_effect[grepl(ra_names, pattern = "\\bz\\b"),1]
g1 <- ra_effect[grepl(ra_names, pattern = "\\bg1\\b"),1]
g2 <- ra_effect[grepl(ra_names, pattern = "\\bg2\\b"),1]

# data wrangling
make_data <- function(thedata,
                      incl.temp = F,
                      incl.cont = F,
                      incl.dens = F) {
  # function to select only complex (nhosts >=3) and with 
  # years with corresponding driver data
  red.data <- dplyr::filter(thedata, nhosts == "3")
  if (incl.temp) red.data <- dplyr::filter(red.data,!is.na(temp_na_rm))
  if (incl.cont) red.data <- dplyr::filter(red.data, !is.na(ord1))
  if (incl.dens) red.data <- dplyr::filter(red.data, !is.na(sfish_density))
  
  
  red.data$slat <- scale(red.data$lat)
  red.data$slat[is.na(red.data$slat)] <- 0
  red.data$slong <- scale(red.data$long)
  red.data$slong[is.na(red.data$slong)] <- 0
  red.data$para.spc <- droplevels(red.data$para.spc)
  red.data$par.type <-droplevels(red.data$par.type)
  red.data$fish.spc <- droplevels(red.data$fish.spc)
  red.data$fish.id <- droplevels(red.data$fish.id)
  red.data$hostpar <- droplevels(red.data$hostpar)
  return(red.data)
}
red.data <- make_data(thedata, incl.temp = T)
red.data <- rem_data(red.data)
U <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$temp_na_rm))
spc <- levels(red.data$para.spc)
nspc <- length(spc)
# create a vector of pararasite species ID for each observation
spcindex <- as.numeric(red.data$para.spc) -1
spc.groups <- levels(red.data$par.type)
n <- as.integer(nrow(red.data))
Xint <- model.matrix(~ -1 + hostpar, data = red.data) 
# separate fish length effects for each parasite / fish combo
Xlength <- Xint * as.numeric(red.data$length)
Xlat <- model.matrix(~ -1 +  para.spc, data = red.data) * as.numeric(red.data$slat) 
Z <- model.matrix(~ -1 + fish.id, data = red.data)

##### Get predicted values #####
logmu <- Xint %*% bo + Xlat %*% blat + Xlength %*% blength + Z %*% z + U %*% g1# + U2 %*% g2# + U3 %*% g3
mu <- exp(logmu)

phi.list <- phi[spcindex]
nsims = 250
mu <- exp(logmu)
##### Simulate data and get scaled residuals #####
simulatedResponse <- matrix(NA, nrow = n, ncol = nsims)
set.seed(123)
for (i in 1:nsims) simulatedResponse[,i] <- rnbinom(n, mu = mu, size = phi.list)
fittedPredictedResponse <- mu
observedResponse <- y

# get scaled residuals
simres <- simulateResidualst(observedResponse, n = 1250, refit = F, integerResponse = T, plot = F, seed = 123, method = c("PIT", "traditional"), rotation = NULL,
                             simulations = simulatedResponse,
                             fitted = fittedPredictedResponse)
##### run DHARMa - like code  #####
scaledResiduals <- simres$scaledResiduals
nObs <- length(observedResponse)
simulationOutput <- list(scaledResiduals = scaledResiduals,
                         observedResponse = observedResponse,
                         fittedPredictedResponse = fittedPredictedResponse,
                         simulatedResponse = simulatedResponse,
                         nObs = length(mu),
                         nSims = nSims,
                         refit = F)

resid.tests_phase2 <- plot.DHARMat(simulationOutput)






