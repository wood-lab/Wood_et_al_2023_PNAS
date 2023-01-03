### Figure 3
### Created by Tim Essington
### essing@uw.edu

## Source code to get expected values
rm(list = ls())
library(TMB)
library(ggplot2)
library(dplyr)
library(tidyr)

load("src/analysis/phase_1_results.Rdata")

compile.and.load <- function(model) {
  compile(paste0("src/", model, ".cpp"))
  return(dyn.load(dynlib(paste0("src/",model))))
}

## Fit four different glmm to parasite data.  Each model has the same base parameters as follows:
## (1) unique intercepts for each parasite / fish species combination
## (2) unique effects of fish size for each parasite / fish species combination
## (3) unique effects of latitude and longitude for each parasite species
## (4) unique effects of temperature and lead for each parasite species
## The models differ with respect to the additional parameters with respect to trends in time
## The are:
## Model 1: no time effects
## Model 3: Species-specific time effects nested within broader taxonomic groupings (2 -level random effect)
## Model 3: Species-specific time effects, nested within taxonomic groupings, with unique taxonomic mean effects based on lifecycle (complex, direct transmission).
## Model 4: Species-specific time effects, nexted within life cycle types

## In all cases, models are fit in two stages.  First a poisson likelihood model is fit because it is more
## numerically stable.  Then, parameters from the poisson likelihood model are used as starting
## values in a negative-binomial likelihood model.
### Setup Data ####

library(tidyr)
library(TMB)
library(ggplot2)

source("src/load_data.R")

# for convenience, make separate design matrix for intercepts (all par / fish combos)
Xint <- model.matrix(~ -1 + hostpar, data = thedata)

# Now make the fixed effects, need to do this in two steps to add in parasite-level fish dependence
# separate fish length effects for each parasite / fish combo
Xlength <- Xint * as.numeric(thedata$length)
Xlat <- model.matrix(~ -1 +  para.spc, data = thedata) * as.numeric(thedata$slat) 
Xlong <- model.matrix(~ -1 + para.spc, data = thedata) * as.numeric(thedata$slong)

# random effects design matrix - parasite specific year effects
U <- model.matrix(~ -1 + para.spc, data = thedata) * as.numeric(scale(thedata$year))
Z <- model.matrix(~ -1 + fish.id, data = thedata)
spc <- levels(thedata$para.spc)
nspc <- length(spc)


# create a vector of parasite species ID for each observation
Xspc <- as.numeric(thedata$para.spc) -1
spc.groups <- levels(thedata$par.type)
y <- thedata$count


# get index of taxonomic group for each species, and create species specific LH type

tax.tmp <- LH.spc.tmp <- NH.spc.tmp <- rep(NA, nspc)
for (i in 1:length(spc)) {
  tax.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$par.type[1]
  LH.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$lifecycle[1]
  NH.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$nhosts[1]
}

# subtract 1 for C++ indexing
taxgroup <- tax.tmp -1
LHspc <- LH.spc.tmp -1
NHspc <- NH.spc.tmp -1

years <- min(thedata$year) : max(thedata$year)
yearsd <- sd(thedata$year)
yearmu = mean(thedata$year)
syears <- (years - yearmu) / yearsd
Xl <- model.matrix(~-1 + nhosts, data = thedata)

### Run TMB model ####
data <- list(n = as.integer(nrow(thedata)),
             nspc = as.integer(nspc),
             y = y,
             U = U,
             Xlat = Xlat,
             Xlength = Xlength,
             Xg = Xl,
             Xint = Xint,
             LH = NHspc,
             Z = Z,
             Xspc = Xspc,
             Upredict = syears,
             npredict = length(syears)
)

parameters <- list(bobar = 0,
                   bo = rep(0, ncol(Xint)),
                   logsigma_bo = 0,
                   bog = rep(0, ncol(Xl)),
                   logsigma_bog = 0,
                   blat = rep(0, ncol(Xlat)),
                   blatbar = 0,
                   logsigma_blat = 0,
                   blength = rep(0, ncol(Xlength)),
                   blengthbar = 0,
                   logsigma_blength = 0,
                   graw = rep(0, ncol(U)),
                   gbar = c(0,0,0),
                   logsigma_g = 0,
                   z = rep(0, ncol(Z)),
                   logsigma_z = 0,
                   logphi = rep(1, nspc),
                   logphibar = 1,
                   loglogsigma_phi = 0
)



model <- "glmm_byspecies_lh_predict"

### Create Handy Function ####
compile.and.load <- function(model) {
  compile(paste0("src/TMB/", model, ".cpp"))
  return(dyn.load(dynlib(paste0("src/TMB/",model))))
}
compile.and.load(model)

obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    random = c("graw", "z", "bo", "bog", "blat", "blength"),
    DLL = model,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)

# extract out the estimated  predicted means and se
est_pars <- summary(rep, "report")
predictions1 <- est_pars[grepl(x = rownames(est_pars), pattern = "\\<mu_predict1\\>"),]
predictions2 <- est_pars[grepl(x = rownames(est_pars), pattern = "\\<mu_predict2\\>"),]
predictions3 <- est_pars[grepl(x = rownames(est_pars), pattern = "\\<mu_predict3\\>"),]

# make a dataframe
all_df <- data.frame(year = years,
                     yhat = c(predictions1[,1],predictions2[,1],predictions3[,1]),
                     se = c(predictions1[,2],predictions2[,2],predictions3[,2]),
                     nhosts = c(rep(1, nrow(predictions1)),
                                rep(2, nrow(predictions2)),
                                rep("3+", nrow(predictions3))
                     )
)
save(all_df,file = "src/analysis/fitted_time_effects.Rdata")
### Make Plot ####
library(ggplot2)
library(dplyr)
load(file = "src/analysis/fitted_time_effects.Rdata")
onehost <- dplyr::filter(all_df, nhosts == 1)
twohost <- dplyr::filter(all_df, nhosts == 2)
threehost <- dplyr::filter(all_df, nhosts == "3+")

library(wesanderson)
cols <- wes_palette(10, name = "Zissou1", type="continuous")[c(7,4,1)]

#### Plot code #####
ggplot(all_df, aes(year, yhat, color = nhosts,
                   ymin = yhat - se, 
                   ymax = yhat + se)) +
  geom_line(size = 2) +
  scale_color_manual(values=rev(cols)) +
  geom_ribbon(data = onehost,
              aes(
                x = year,
                ymin = yhat - se,
                ymax = yhat + se
              ),
              inherit.aes = FALSE,
              fill =  cols[3],
              alpha = 0.4,
              show.legend = F
  ) +
  geom_ribbon(
    data = twohost,
    aes(
      x = year,
      ymin = yhat - se,
      ymax = yhat + se
    ),
    inherit.aes = FALSE,
    fill = cols[2],
    alpha = 0.4,
    show.legend = F
  ) +
  geom_ribbon(
    data = threehost,
    aes(
      x = year,
      ymin = yhat - se,
      ymax = yhat + se
    ),
    inherit.aes = FALSE,
    fill = cols[1],
    alpha = 0.4,
    show.legend = F
  ) +
  labs(color = "") +
  #  geom_ribbon(alpha = 0.4, outline.type ="both", show.legend = F) + 
  scale_y_continuous(limits = c(0, 1.65), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "year", y = "parasite count") + 
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 13, color = "black")) +
  theme(axis.title= element_text(size = 20)) +
  theme(legend.text = element_text(size = 30, color = "black"))+
  theme(legend.position = c(0.7,0.75))

# save image in figures
ggsave("figures/fitted_time_effects.jpeg", width = 5, height = 4, units = "in")

