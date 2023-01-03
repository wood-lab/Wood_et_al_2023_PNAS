### Phase 2 analysis
### Created by Tim Essington
### essing@uw.edu

## Fit 4 different glmm to parasite data, only including those with complex life history.
#.  Each model has the same base parameters as follows:
## (1) unique intercepts for each parasite / fish species combination
## (2) unique effects of fish size for each parasite / fish species combination
## (3) unique effects of latitude and longitude for each parasite species
## (4) random effect for fish specimen


## Then consider as additional parameters
## (A) Species-specific time effects (RE)
## (B) Species-specific contaminant effects (RE)
## (C) Species-specific temperature effects (RE)
## (D) Species-specific fish density effect (RE)

## In all cases, models are fit in two stages.  First a poisson likelihood model is fit because it is more
## numerically stable.  Then, parameters from the poisson likelihood model are used as starting
## values in a negative-binomial likelihood model.

### Setup Data ####

library(tidyr)
library(dplyr)
library(TMB)
library(ggplot2)
library(viridis)

source("src/load_data.R")

# get phase 1 estimates
load("src/analysis/phase_1_results.Rdata")
# get SE and MLE of 
rep.2.use <- models$rep.byspecieslh
fixed.year.effect <- summary(rep.2.use, "fixed")
year.effect <- fixed.year.effect[grep(rownames(fixed.year.effect), pattern = "gbar"),][3,]


## data filter - remove any NAs, also only pull out the complex life cycle
# function to remove data if there are not sufficient positive occurences

### Setup Functions ####
rem_data <- function(red.data) {
  notzero <- function(x)
    min(1, x)
  tmp.data <- red.data
  red.data$foo <- sapply(X = red.data$count, FUN = notzero)
  data_summary <- red.data %>%
    group_by(para.spc, fish.spc) %>%
    summarise(ntotal = n(), nfoo = sum(foo))
  
  data_summary$prop <- data_summary$nfoo / data_summary$ntotal
  if (min(data_summary$prop) < 0.05) {
    data_rem <- data_summary[data_summary$prop < 0.05, ]
    rem.par.spc <- as.array(as.character(data_rem$para.spc))
    rem.fish.spc <- as.array(as.character(data_rem$fish.spc))
    
    tmp.data$para.spc <- as.character(tmp.data$para.spc)
    tmp.data$fish.spc <- as.character(tmp.data$fish.spc)
    for (i in 1:nrow(data_rem)) {
      index <-
        which(tmp.data$para.spc == rem.par.spc[i] &
                tmp.data$fish.spc == rem.fish.spc[i])
      tmp.data <- tmp.data[-index, ]
    }
    tmp.data$para.spc <- as.factor(tmp.data$para.spc)
    tmp.data$fish.spc <- as.factor(tmp.data$fish.spc)
  }
  return(tmp.data)
}

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


# get the SD of year, used to convert slopes to per year
year.scale <- sd(thedata$year)

compile.and.load <- function(model) {
  compile(paste0("src/TMB/", model, ".cpp"))
  return(dyn.load(dynlib(paste0("src/TMB/",model))))
}


fit_model <- function(U,red.data, npred = 1, U2 = NULL, U3 = NULL) {
  spc <- levels(red.data$para.spc)
  nspc <- length(spc)
  # create a vector of pararasite species ID for each observation
  spcindex <- as.numeric(red.data$para.spc) -1
  spc.groups <- levels(red.data$par.type)
  y <- red.data$count
  n <- as.integer(nrow(red.data))
  Xint <- model.matrix(~ -1 + hostpar, data = red.data) 
  # separate fish length effects for each parasite / fish combo
  Xlength <- Xint * as.numeric(red.data$length)
  Xlat <- model.matrix(~ -1 +  para.spc, data = red.data) * as.numeric(red.data$slat) 
  Z <- model.matrix(~ -1 + fish.id, data = red.data)
  
  ##### Fit model with no predictor ####
  if (npred ==0) {
    data <- list(n = as.integer(nrow(red.data)),
                 nspc = as.integer(nspc),
                 y = y,
                 Xlat = Xlat,
                 Xlength = Xlength,
                 Xint = Xint,
                 Z = Z,
                 spcindex = as.integer(spcindex))
    
    
    parameters <- list(bo = rep(0,ncol(Xint)),
                       logsigma_bo = 0,
                       bobar = 0,
                       blat = rep(0, ncol(Xlat)),
                       blatbar = 0,
                       logsigma_blat = 0,
                       blength = rep(0, ncol(Xlength)),
                       blengthbar = 0,
                       logsigma_blength = 0,
                       z = rep(0, ncol(Z)),
                       logsigma_z = 0,
                       logphi = rep(1, nspc),
                       logphibar = 1,
                       loglogsigma_phi = 0
    )
    
    
    model <- "glmm_no_pred"
    compile.and.load(model)
    obj <-
      MakeADFun(
        data = data,
        parameters = parameters,
        DLL = model,
        random = c("z", "bo", "blength", "blat", "logphi"),
        silent = TRUE
      )
  }
  ##### Fit model with one predictor ####
  
  if (npred ==1) {
    data <- list(n = as.integer(nrow(red.data)),
                 nspc = as.integer(nspc),
                 y = y,
                 U =U,
                 Xlat = Xlat,
                 Xlength = Xlength,
                 Xint = Xint,
                 Z = Z,
                 spcindex = as.integer(spcindex))
    
    
    parameters <- list(bo = rep(0,ncol(Xint)),
                       logsigma_bo = 0,
                       bobar = 0,
                       blat = rep(0, ncol(Xlat)),
                       blatbar = 0,
                       logsigma_blat = 0,
                       blength = rep(0, ncol(Xlength)),
                       blengthbar = 0,
                       logsigma_blength = 0,
                       g = rnorm(ncol(U), 0, 1),
                       gbar = 0,
                       logsigma_g = 0,
                       z = rep(0, ncol(Z)),
                       logsigma_z = 0,
                       logphi = rep(1, nspc),
                       logphibar = 1,
                       loglogsigma_phi = 0
    )
    
    
    model <- "glmm_byspecies"
    compile.and.load(model)
    obj <-
      MakeADFun(
        data = data,
        parameters = parameters,
        DLL = model,
        random = c("z", "bo", "g", "blength", "blat", "logphi"),
        silent = TRUE
      )
  }
  ##### Fit model with two predictors ####
  if (npred == 2) {
    data <- list(n = n,
                 nspc = as.integer(nspc),
                 y = y,
                 U1 = U,
                 U2 = U2,
                 Xlat = Xlat,
                 Xlength = Xlength,
                 Xint = Xint,
                 spcindex = as.integer(spcindex),
                 Z = Z
    )
    
    parameters <- list(bo = rep(0,ncol(Xint)),
                       bobar = 0,
                       logsigma_bo = 0,
                       blat = rep(0, ncol(Xlat)),
                       blatbar = 0,
                       logsigma_blat = 0,
                       blength = rep(0, ncol(Xlength)),
                       blengthbar = 0,
                       logsigma_blength = 0,
                       g1 = rep(0, ncol(U)),
                       gbar1 = 0,
                       logsigma_g1 = 0,
                       g2 = rep(0, ncol(U2)),
                       gbar2 = 0,
                       logsigma_g2 = 0,
                       z = rep(0, ncol(Z)),
                       logsigma_z = 0,
                       logphi = rep(1, nspc),
                       logphibar = 1,
                       loglogsigma_phi = 0
    )
    
    model <- "glmm_byspecies_two_predict"
    compile.and.load(model)
    obj <-
      MakeADFun(
        data = data,
        parameters = parameters,
        random = c("z", "bo", "blat", "blength", "g1", "g2", "logphi"),
        DLL = "glmm_byspecies_two_predict",
        silent = TRUE,
        checkParameterOrder = TRUE
      )
    
  }
  #### Fit model with three predictors ####
  if (npred == 3) {
    data <- list(n = n,
                 nspc = as.integer(nspc),
                 y = y,
                 U1 = U,
                 U2 = U2,
                 U3 = U3,
                 Xlat = Xlat,
                 Xlength = Xlength,
                 Xint = Xint,
                 spcindex = as.integer(spcindex),
                 Z = Z
    )
    
    parameters <- list(bo = rep(0,ncol(Xint)),
                       bobar = 0,
                       logsigma_bo = 0,
                       blat = rep(0, ncol(Xlat)),
                       blatbar = 0,
                       logsigma_blat = 0,
                       blength = rep(0, ncol(Xlength)),
                       blengthbar = 0,
                       logsigma_blength = 0,
                       g1 = rep(0, ncol(U)),
                       gbar1 = 0,
                       logsigma_g1 = 0,
                       g2 = rep(0, ncol(U2)),
                       gbar2 = 0,
                       logsigma_g2 = 0,
                       g3 = rep(0, ncol(U3)),
                       gbar3 = 0,
                       logsigma_g3 = 0,
                       z = rep(0, ncol(Z)),
                       logsigma_z = 0,
                       logphi = rep(1, nspc),
                       logphibar = 1,
                       loglogsigma_phi = 0
    )
    
    model <- "glmm_byspecies_three_predict"
    compile.and.load(model)
    obj <-
      MakeADFun(
        data = data,
        parameters = parameters,
        random = c("z", "bo", "blat", "blength", "g1", "g2","g3", "logphi"),
        DLL = "glmm_byspecies_three_predict",
        silent = TRUE,
        checkParameterOrder = TRUE
      )
    
  }
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep = sdreport( obj,
                  getReportCovariance = FALSE
  )
  
  return(list(opt = opt, rep = rep, obj = obj, data = data))
}

TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}

models <- list()

### Run Model 1: Temperature ####
red.data <- make_data(thedata, incl.temp = T)
red.data <- rem_data(red.data)

null.model <- fit_model(U = NULL, npred = 0, red.data = red.data)
models$temp$null <- null.model

U <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$temp_na_rm))
temp.model <- fit_model(U = U, red.data = red.data)
models$temp$temp <- temp.model
Uyear <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$year))
year.model <- fit_model(U = Uyear, red.data = red.data)
models$temp$year <- year.model
U2 <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$year))
U1 <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$temp_na_rm))
temp.year.model <- fit_model(U = U1, U2 = U2, npred = 2, red.data = red.data)
models$temp$tempyear <- temp.year.model


re <- summary(temp.year.model$rep, "random")
g1 <- re[grep(rownames(re), pattern = "\\bg1\\b"),]
g2 <- re[grep(rownames(re), pattern = "\\bg2\\b"),]
## Prepare plot for analysis
gs <- tibble(g.t.hat = g1[,1],
             g.t.se = g1[,2],
             g.y.hat = g2[,1],
             g.y.se = g2[,2]
)

col <- viridis(n = 10, option = "mako") [3]
colse <- viridis(n = 10, option = "mako") [7]
ggplot(gs, aes(x = g.t.hat, y = g.y.hat)) +
  geom_errorbar(aes(ymin = g.y.hat - g.y.se,
                    ymax = g.y.hat + g.y.se),
                size = 0.35,
                color = colse) + 
  geom_errorbar(aes(xmin = g.t.hat - g.t.se,
                    xmax = g.t.hat + g.t.se),
                size = 0.35,
                color = colse) + 
  geom_point(colour = col, size = 3.5) + 
  # ylim(-1.7, 1.7) + 
  #  xlim(-0.5, 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab(bquote(gamma[k*","*temp])) +
  ylab(bquote(gamma[k*","*year])) +
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())

#### Run model with smoothed temperature
U <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$temp_smooth))
smooth.temp.model <- fit_model(U = U, red.data = red.data)

##### create AIC table for Temperature vs. Year Comparison ####
AIC <- rep(NA, 4)
AIC[1] <- TMBAIC(null.model$opt, p = 2, n = nrow(red.data))
AIC[2] <- TMBAIC(temp.model$opt, p = 2, n = nrow(red.data))
AIC[3] <- TMBAIC(year.model$opt, p = 2, n = nrow(red.data))
AIC[4] <- TMBAIC(temp.year.model$opt, p = 2, n = nrow(red.data))
#AIC[5] <- TMBAIC(smooth.temp.model$opt, p = 2, n = nrow(red.data))

dAIC <- AIC - min(AIC)
names(dAIC) <- c("null","temp", "year", "temp + year")
print(dAIC)





### Run Model 2 : Both Contaminant ordinates ####
red.data <- make_data(thedata, incl.cont = T)
red.data <- rem_data(red.data)

null.model <- fit_model(U = NULL, npred = 0, red.data = red.data)
models$cont$null <- null.model
U1 <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$ord1raw))
U2 <-model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$ord2raw))
cont.model <- fit_model(U = U1, red.data = red.data, npred = 2, U2 = U2)
models$cont$cont <- cont.model
obj <- cont.model$obj
print(paste0("contaminant model nll = ", cont.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_cont <- r$nll_data

Uyear <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$year))
year.model <- fit_model(U = Uyear, red.data = red.data)
models$cont$year <- year.model

obj <- year.model$obj
print(paste0("year model nll = ", year.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_year <- r$nll_data


year.cont.model <- fit_model(U = U1, red.data = red.data, npred = 3, U2 = U2, U3 = Uyear)
models$cont$contyear <- year.cont.model
obj <- year.cont.model$obj
print(paste0("year + contaminant model nll = ", year.cont.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_year_cont <- r$nll_data

re <- summary(year.cont.model$rep, "random")
g1 <- re[grep(rownames(re), pattern = "\\bg1\\b"),]
g2 <- re[grep(rownames(re), pattern = "\\bg3\\b"),]
## Prepare plot for analysis
gs <- tibble(g.t.hat = g1[,1],
             g.t.se = g1[,2],
             g.y.hat = g2[,1],
             g.y.se = g2[,2]
)
library(viridis)
library(ggplot2)
col <- viridis(n = 10, option = "mako") [3]
colse <- viridis(n = 10, option = "mako") [7]
ggplot(gs, aes(x = g.t.hat, y = g.y.hat)) +
  geom_errorbar(aes(ymin = g.y.hat - g.y.se,
                    ymax = g.y.hat + g.y.se),
                size = 0.35,
                color = colse) + 
  geom_errorbar(aes(xmin = g.t.hat - g.t.se,
                    xmax = g.t.hat + g.t.se),
                size = 0.35,
                color = colse) + 
  geom_point(colour = col, size = 3.5) + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab(bquote(gamma[k*","*cont_ord1])) +
  ylab(bquote(gamma[k*","*year])) +
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())


##### create AIC table for Contaminant vs. Year Comparison ####
AIC <- rep(NA, 4)
AIC[1] <- TMBAIC(null.model$opt, n = nrow(red.data))
AIC[2] <- TMBAIC(cont.model$opt, n = nrow(red.data))
AIC[3] <- TMBAIC(year.model$opt, n = nrow(red.data))
AIC[4] <- TMBAIC(year.cont.model$opt, n = nrow(red.data))

dAIC <- AIC - min(AIC)
names(dAIC) <- c("Null","Cont", "Year", "Cont + Year")
print(dAIC)


### Run Model 3: fish density ####
red.data <- make_data(thedata, incl.dens = T)
red.data <- rem_data(red.data)
null.model <- fit_model(U = NULL, npred = 0, red.data = red.data)
models$density$null <- null.model
obj <- null.model$obj
print(paste0("null model nll = ", null.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_null <- r$nll_data


U <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(red.data$sfish_density)
density.model <- fit_model(U = U, red.data = red.data)
models$density$density <- density.model
obj <- density.model$obj
print(paste0("density model nll = ", density.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_density<- r$nll_data

Uyear <- model.matrix(~ -1 + para.spc, data = red.data) * as.numeric(scale(red.data$year))
year.model <- fit_model(U = Uyear, npred = 1, red.data = red.data)
models$density$year <- year.model
obj <- year.model$obj
print(paste0("year model nll = ", year.model$opt$objective))
lp <- obj$env$last.par.best
r <- obj$report(lp)
nll_data_year <- r$nll_data
year.density.model <- fit_model(U= U, npred = 2, red.data = red.data, U2 = Uyear)
models$density$densityyear <- year.density.model


re <- summary(density.model$rep, "random")
g1 <- re[grep(rownames(re), pattern = "\\bg\\b"),]
re <- summary(year.model$rep, "random")
g2 <- re[grep(rownames(re), pattern = "\\bg\\b"),]
## Prepare plot for analysis
gs <- tibble(g.t.hat = g1[,1],
             g.t.se = g1[,2],
             g.y.hat = g2[,1],
             g.y.se = g2[,2]
)

col <- viridis(n = 10, option = "mako") [3]
colse <- viridis(n = 10, option = "mako") [7]
ggplot(gs, aes(x = g.t.hat, y = g.y.hat)) +
  geom_errorbar(aes(ymin = g.y.hat - g.y.se,
                    ymax = g.y.hat + g.y.se),
                size = 0.35,
                color = colse) + 
  geom_errorbar(aes(xmin = g.t.hat - g.t.se,
                    xmax = g.t.hat + g.t.se),
                size = 0.35,
                color = colse) + 
  geom_point(colour = col, size = 3.5) + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab(bquote(gamma[k*","*density])) +
  ylab(bquote(gamma[k*","*year])) +
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())

AIC <- rep(NA, 4)
AIC[1] <- TMBAIC(null.model$opt, n = nrow(red.data))
AIC[2] <- TMBAIC(density.model$opt, n = nrow(red.data))
AIC[3] <- TMBAIC(year.model$opt, n = nrow(red.data))
AIC[4] <- TMBAIC(year.density.model$opt, n = nrow(red.data))

dAIC <- AIC - min(AIC)
names(dAIC) <- c("Null", "Density", "Year", "Density + Year")
print(dAIC)

# ### Summarize and Save output ####
get_gbar <- function(rep) {
  fixed.par <- summary(rep, "fixed")
  gbar <- fixed.par[grepl(rownames(fixed.par), pattern = "gbar"),1]
  gbarse <- fixed.par[grepl(rownames(fixed.par), pattern = "gbar"),2]
  return(cbind(gbar, gbarse))
}
#   
#   
# gbar.year <- get_gbar(year.model$rep)
gbar.temp <- get_gbar(temp.model$rep)
gbar.cont <- get_gbar(cont.model$rep)
gbar.density <- get_gbar(density.model$rep)

fit_summary <- tibble(predictor  = c("year","temperature", "contaminant PC1", "contaminant PC2", "density"),
                      slope = c(year.effect[1],gbar.temp[1], gbar.cont[,1], gbar.density[1]),
                      se = c(year.effect[2],gbar.temp[2], gbar.cont[,2], gbar.density[2])
)
# 
# 
# 
models$fit_summary <- fit_summary
save(models, file = "src/analysis/phase_2_results.Rdata")
# 
# load(file = "src/analysis/phase_2_results.Rdata")
# # ### make data tables ####
#  for (i in 1:4) {
#    
#      model.2.use <- models[[i]]
#      fixedef <- summary(model.2.use$rep, "fixed")
#      report <- summary(model.2.use$rep, "report") # for getting sigmas
#      
#      fixed_names <- rownames(fixedef)
#      rep_names <- rownames(report)
#      if (!i==3) {
#       coef.report <- data.frame(parameter = c("gbar", "sigma_g", "bobar","sigma_bo", "blatbar", "sigma_blat",
#                                               "blengthbar", "sigma_blength","sigma_z"),
#                                Est = NA,
#                                SE = NA)
#       # extract means
#       for (j in c(1,3,5,6,7)) {
#         coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),1]
#         coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),2]
#       }
#       
#       # extract sigmas
#       for (j in c(2,4,6,8,9)) {
#         coef.report$Est[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),1]
#         coef.report$SE[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),2]
#       }
#       print(paste0("model ", i))
#       print(coef.report)
#      }
#      if (i ==3) {
#        coef.report <- data.frame(parameter = c("gbar1", "sigma_g1", "gbar2", "sigma_g2", "bobar","sigma_bo", "blatbar", "sigma_blat",
#                                                "blengthbar", "sigma_blength","sigma_z"),
#                                  Est = NA,
#                                  SE = NA)
#        # extract means
#        for (j in c(1,3,5,6,7,9)) {
#          coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),1]
#          coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),2]
#        }
#        
#        # extract sigmas
#        for (j in c(2,4,6,8,10,11)) {
#          coef.report$Est[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),1]
#          coef.report$SE[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),2]
#        }
#        print(paste0("model ", i))
#           print(coef.report)
#      }
#  }
### Create Parameter Estimation Table ####
get_gbar <- function(mod, pattern = "gbar") {
  fe <- summary(mod$rep, "fixed")
  gbar <-unname(fe[grep(rownames(fe), pattern = pattern),])
  return(gbar)
}

make_coefs <- function(model, is.cont = F) {
  model_names <- names(model)
  if(is.cont) {
    coefs <- data_frame(model = c("pred", "year","pred+year"),
                        pred1 = NA,
                        pred1se = NA,
                        pred2 = NA,
                        pred2se = NA,
                        year = NA,
                        yearse = NA)
    
    gbar <- get_gbar(model[[2]])
    gbars <- matrix(c(gbar[1,1:2], gbar[2,1:2]), nrow =1, ncol = 4)
    coefs[1,2:5] <- gbars
    gbar <- get_gbar(model[[3]])
    coefs[coefs$model == "year",6:7] <- matrix(gbar, nrow =1, ncol = 2)
    gbar <- get_gbar(model[[4]])
    gbars <- matrix(c(gbar[1,1:2],gbar[3,1:2], gbar[2,1:2]), nrow =1, ncol = 6)
    coefs[coefs$model == "pred+year", 2:7] <- gbars
    
  }
  if(!is.cont) {
    coefs <- data_frame(model = c("pred", "year","pred+year"),
                        pred = NA,
                        predse = NA,
                        year = NA,
                        yearse = NA)
    gbar <- get_gbar(model[[2]])
    coefs[1,2] <- gbar[1]
    coefs[1,3] <- gbar[2]
    gbar <- get_gbar(model[[3]])
    coefs[2,4] <- gbar[1]
    coefs[2,5] <- gbar[2]
    gbar <- get_gbar(model[[4]])
    gbars <- matrix(c(gbar[1,1:2], gbar[2,1:2]), nrow =1, ncol = 4)
    coefs[coefs$model == "pred+year", 2:5] <- gbars
  }
  return(coefs)
}
#### Temperature ####
make_coefs(models$temp)
make_coefs(models$cont, is.cont= T)
make_coefs(models$density)


