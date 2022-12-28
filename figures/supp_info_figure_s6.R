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
load("src/analysis/phase_1_results.Rdata")

# extract parameters
fixed_effect <- summary(models$rep.byspecieslh, "fixed")
fixed_names <- rownames(fixed_effect)
rep_effect <- summary(models$rep.byspecieslh, "report")
rep_names <- rownames(rep_effect)
ra_effect <- summary(models$rep.byspecieslh, "random")
ra_names <- rownames(ra_effect)


logphi <- ra_effect[grepl(ra_names, pattern = "logphi"),1]
phi <- exp(logphi)
bo <-  ra_effect[grepl(ra_names, pattern = "\\bbo\\b"),1]
bog <-  ra_effect[grepl(ra_names, pattern = "\\bbog\\b"),1]
blat<-  ra_effect[grepl(ra_names, pattern = "\\bblat\\b"),1]
blength<-  ra_effect[grepl(ra_names, pattern = "\\bblength\\b"),1]
gbar <- fixed_effect[grepl(fixed_names, pattern = "\\bgbar\\b"),1]
z <- ra_effect[grepl(ra_names, pattern = "\\bz\\b"),1]
spc_slopes <- rep_effect[grepl(rep_names, pattern = "spc_slope"),1]
# create a vector of parasite species ID for each observation
Xspc <- as.numeric(thedata$para.spc)
spc.groups <- levels(thedata$par.type)
y <- thedata$count

##### Get predicted values #####
logmu <- Xint %*% bo + Xg %*% bog + Xlat %*% blat + Xlength %*% blength + Z %*% z + U %*% spc_slopes
mu <- exp(logmu)
phi.list <- phi[Xspc]
nsims = 250

##### Simulate data and get scaled residuals #####
simulatedResponse <- matrix(NA, nrow = nrow(thedata), ncol = nsims)
set.seed(123)
for (i in 1:nSims) simulatedResponse[,i] <- rnbinom(nrow(thedata), mu = mu, size = phi.list)
fittedPredictedResponse <- mu
observedResponse <- thedata$count

# get scaled residuals
simres <- simulateResidualst(observedResponse, n = 250, refit = F, integerResponse = T, plot = F, seed = 123, method = c("PIT", "traditional"), rotation = NULL,
                             simulations = simulatedResponse,
                             fitted = fittedPredictedResponse)
##### run DHARMa - like code  #####
scaledResiduals <- simres$scaledResiduals

simulationOutput <- list(scaledResiduals = scaledResiduals,
                         observedResponse = observedResponse,
                         fittedPredictedResponse = fittedPredictedResponse,
                         simulatedResponse = simulatedResponse,
                         nObs = length(mu),
                         nSims = nSims,
                         refit = F)

resid.tests <- plot.DHARMat(simulationOutput)
