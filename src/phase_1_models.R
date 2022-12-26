library(tidyr)
library(TMB)
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)

### Setup Data ####
source("src/load_data.R")
source("src/configure_TMB_runs.R")                       



### Create Handy Function ####
compile.and.load <- function(model) {
  compile(paste0("src/TMB/", model, ".cpp"))
  return(dyn.load(dynlib(paste0("src/TMB/",model))))
}
### Run Model 1 ####
data <- list(n = as.integer(nrow(thedata)),
             nspc = as.integer(nspc),
             y = y,
             Xlat = Xlat,
             Xlength = Xlength,
             Xint = Xint,
             Z = Z,
             Xspc = as.integer(Xspc))


parameters <- list(bo = rep(0,ncol(Xint)),
                   bobar = 0,
                   logsigma_bo = 0,
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


model <- "glmm"
compile.and.load(model)
obj <- MakeADFun(data = data, parameters = parameters, DLL = model, random = c("z","bo","blat", "blength","logphi"), silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
opt.null <- opt
rep.null <- rep

### Run Model 2 ####

data <- list(n = as.integer(nrow(thedata)),
             nspc = as.integer(nspc),
             y = y,
             U =U,
             Xlat = Xlat,
             Xlength = Xlength,
             Xint = Xint,
             Z = Z,
             spcindex = as.integer(Xspc))


parameters <- list(bo = rep(0,ncol(Xint)),
                   bobar = 0,
                   logsigma_bo = 0,
                   blat = rep(0, ncol(Xlat)),
                   blatbar = 0,
                   logsigma_blat = 0,
                   blength = rep(0, ncol(Xlength)),
                   blengthbar = 0,
                   logsigma_blength = 0,
                   g = rep(0, ncol(U)),
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
    random = c("z", "bo", "blat", "blength", "g", "logphi"),
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
opt.byspecies <- opt
rep.byspecies <- rep

## Run Model 3 ####
# This uses the number of hosts as categories, instead of direct vs. complex. 

data <- list(n = as.integer(nrow(thedata)),
             nspc = as.integer(nspc),
             y = y,
             U = U,
             Xlat = Xlat,
             Xlength = Xlength,
             Xint = Xint,
             Xg = Xg,
             Xspc = Xspc,
             LH = NHspc,
             Z = Z)

parameters <- list(bo = rep(0,ncol(Xint)),
                   bobar = 0,
                   logsigma_bo = 0,
                   bog = rep(0, ncol(Xg)),
                   logsigma_bog = 0,
                   blat = rep(0, ncol(Xlat)),
                   blatbar = 0,
                   logsigma_blat = 0,
                   blength = rep(0, ncol(Xlength)),
                   blengthbar = 0,
                   logsigma_blength = 0,
                   graw = rep(0, ncol(U)),
                   gbar =rep(0, 3),
                   logsigma_g = 0,
                   z = rep(0, ncol(Z)),
                   logsigma_z = 0,
                   logphi = rep(1,nspc),
                   logphibar = 1,
                   loglogsigma_phi = 0
)


model <- "glmm_byspecies_lh"
compile.and.load(model)

obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    random = c("z","bo","bog","blat", "blength", "graw","logphi"),
    DLL = model,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
opt.byspecieslh <- opt
rep.byspecieslh <- rep


### Run Model 4 ####
data <- list(n = as.integer(nrow(thedata)),
             nspc = as.integer(nspc),
             y = y,
             U = U,
             Xlat = Xlat,
             Xlength = Xlength,
             Xint = Xint,
             Xg = Xa,
             Xspc = as.integer(Xspc),
             A = A,
             Z = Z
)

parameters <-
  list(
    bo = rep(0,ncol(Xint)),
    bobar = 0,
    logsigma_bo = 0,
    bog = rep(0, ncol(Xa)),
    logsigma_bog = 0,
    blat = rep(0, ncol(Xlat)),
    blatbar = 0,
    logsigma_blat = 0,
    blength = rep(0, ncol(Xlength)),
    blengthbar = 0,
    logsigma_blength = 0,
    g = rep(0, ncol(U)),
    a = rep(0, ncol(A)),
    logsigma_g = 0,
    z = rep(0, ncol(Z)),
    logsigma_z = 1,
    logphi = rep(1,nspc),
    logphibar = 1,
    loglogsigma_phi = 0
  )

model <- "glmm_byspecies_groups"
compile.and.load(model)
obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    random = c("z","bo", "bog","blat", "blength", "g","logphi"),
    DLL = "glmm_byspecies_groups",
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = FALSE)
opt.byspeciesgroups <- opt
rep.byspeciesgroups <- rep

### Save Model Results ####
filename <- "src/analysis/phase_1_results.Rdata"
models = list(opt.null = opt.null,
              rep.null = rep.null,
              opt.byspecies = opt.byspecies,
              rep.byspecies = rep.byspecies,
              opt.byspecieslh = opt.byspecieslh,
              rep.byspecieslh = rep.byspecieslh,
              opt.byspeciesgroups = opt.byspeciesgroups,
              rep.byspeciesgroups = rep.byspeciesgroups)
save(models, file = filename)


### Model Selection ####
TMBAIC=function(opt, p=2, n=Inf){
  k = length(opt[["par"]])
  if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
  if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
  Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
  return( Return )
}
load(file = "src/analysis/phase_1_results.Rdata")
AIC <- rep(NA, 4)
AIC[1] = TMBAIC(models$opt.null, p = 2, n = nrow(thedata))
AIC[2] = TMBAIC(models$opt.byspecies, p = 2, n = nrow(thedata))
AIC[3] = TMBAIC(models$opt.byspecieslh, p = 2, n = nrow(thedata))
AIC[4] = TMBAIC(models$opt.byspeciesgroups, p = 2, n = nrow(thedata))

DAIC <- AIC - min(AIC)
print(DAIC)

AIC_table <- tibble(model =c ( "null", "time", "time by number of hosts", "time by phyla"),
                    npars = c(length(models$opt.null$par),
                              length(models$opt.byspecies$par),
                              length(models$opt.byspecieslh$par),
                              length(models$opt.byspeciesgroups$par)
                    ), 
                    NLL = c(models$opt.null$objective,
                            models$opt.byspecies$objective,
                            models$opt.byspecieslh$objective,
                            models$opt.byspeciesgroups$objective
                    ),
                    DAIC = num(DAIC,
                               digits = 2)
)


print(AIC_table)


### Plot random effects histograms ####
# no grouping
rep.byspecies <- models$rep.byspecies
randEffects <- summary(rep.byspecies, "report")
randEffectsNames <- rownames(randEffects)
gs.byspecies <- randEffects[grepl(randEffectsNames, pattern = "spc_slopes"),1]
tau.byspecies<- 1/randEffects[grepl(randEffectsNames, pattern = "spc_slopes"),2]


# grouping
rep.byspecieslh <- models$rep.byspecieslh
randEffects <- summary(rep.byspecieslh, "report")
randEffectsNames <- rownames(randEffects)
gs.byspeciesgroups <- randEffects[grepl(randEffectsNames, pattern = "spc_slopes"),1]
tau.byspeciesgroups <- 1/randEffects[grepl(randEffectsNames, pattern = "spc_slopes"),2]

plotdf <-
  data.frame(model = c(
    rep(
      "No grouping",
      times = length(gs.byspecies)),
    rep("Group by # hosts", times = length(gs.byspecies))
  ),
  ghat = c(gs.byspecies, gs.byspeciesgroups),
  tau = c(tau.byspecies, tau.byspeciesgroups),
  nhost = rep(c("1", "2", "3+")[NHspc + 1], times = 2)
  )

plotdf$model  <- as.factor(plotdf$model)

plotdf %>%
  mutate(model = fct_relevel(model, "No grouping", "Group by # hosts"))

png(filename = "src/figures/compare_ranef.png",
    width = 6,
    height = 4,
    units = "in",
    res = 300
)
plotdf %>%
  mutate(model = fct_relevel(model, "No grouping", "Group by # hosts")) %>% 
  ggplot(aes(x = ghat, y = ..density.., weight = tau)) +
  geom_histogram(bins = 15, colour = "black", fill = "gray") +
  facet_grid(rows = vars(model), cols = vars(nhost), scales = "free", labeller = label_context) +
  labs(x = "Random Effect", y = "Density") +
  scale_x_continuous( limits = c(-1.25, 1.25), expand = c(0,0)) + 
  scale_y_continuous( limits = c(0, 1.55), expand = c(0,0)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = "black")) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title= element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) + 
  theme(strip.background = element_blank()) + 
  theme(axis.line=element_line())
dev.off()                               


### Make parameter tables ####
for (i in 1:4) {
  model.2.use <- models[[2*i]]
  fixedef <- summary(model.2.use, "fixed")
  report <- summary(model.2.use, "report") # for getting sigmas
  
  fixed_names <- rownames(fixedef)
  rep_names <- rownames(report)
  ###### For Null #####
  if (i ==1) { 
    coef.report <- data.frame(parameter = c("bobar","sigma_bo", "blatbar", "sigma_blat",
                                            "blengthbar", "sigma_blength","sigma_z",
                                            "logphibar", "logsigma_phi"),
                              Est = NA,
                              SE = NA)
    
    for (j in c(1,3,5,8)) {
      coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),1]
      coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),2]
    }
    
    # extract sigmas
    for (j in c(2,4,6,7,9)) {
      coef.report$Est[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),1]
      coef.report$SE[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),2]
    }
    print(paste0("model ", i))
    print(coef.report)  
  }
  ### For Year ####
  if (i ==2) { 
    coef.report <- data.frame(parameter = c("gbar", "sigma_g","bobar","sigma_bo", "blatbar", "sigma_blat",
                                            "blengthbar", "sigma_blength","sigma_z",
                                            "logphibar", "logsigma_phi"),
                              Est = NA,
                              SE = NA)
    for (j in c(1,3,5,7,10)) {
      coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),1]
      coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = coef.report[j,1]),2]
    }
    
    # extract sigmas
    for (j in c(2,4,6,8,9,11)) {
      coef.report$Est[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),1]
      coef.report$SE[j] <- report[grepl(rep_names, pattern = coef.report[j,1]),2]
    }
    print(paste0("model ", i))
    print(coef.report)  
  }
  
  
  ### For Year by lifecycle ####
  if (i ==3) {
    coef.report <- data.frame(parameter = c("gbar1", "gbar2", "gbar3", "sigma_g", "bobar","sigma_bo", "sigma_bog", "blatbar", "sigma_blat",
                                            "blengthbar", "sigma_blength","sigma_z","logphibar", "logsigma_phi"),
                              Est = NA,
                              SE = NA)
    # extract means for gpar
    coef.report$Est[c(1,2,3)] <- fixedef[grepl(fixed_names, pattern = "gbar"),1]
    coef.report$SE[c(1,2,3)] <- fixedef[grepl(fixed_names, pattern = "gbar"),2]
    
    # extract means
    for (j in c(5,8,10,13)) {
      parname <- paste0("\\b",coef.report[j,1], "\\b")
      coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = parname),1]
      coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = parname),2]
    }
    
    # extract sigmas
    for (j in c(4,6,7,9,11,12,14)) {
      parname <- paste0("\\b",coef.report[j,1], "\\b")
      coef.report$Est[j] <- report[grepl(rep_names, pattern = parname),1]
      coef.report$SE[j] <- report[grepl(rep_names, pattern = parname),2]
    }
    print(paste0("model ", i))
    print(coef.report)
  }
  ### For Year by taxa ####
  if (i ==4) {
    coef.report <- data.frame(parameter = c(rep("a",7), "sigma_g", "bobar","sigma_bo", "sigma_bog", "blatbar", "sigma_blat",
                                            "blengthbar", "sigma_blength","sigma_z","logphibar", "logsigma_phi"),
                              Est = NA,
                              SE = NA)
    
    # extract means for gpar
    coef.report$Est[1:7] <- fixedef[grepl(fixed_names, pattern = "\\ba\\b"),1]
    coef.report$SE[1:7] <- fixedef[grepl(fixed_names, pattern = "\\ba\\b"),2]
    # extract means
    for (j in c(9,12,14,17)) {
      parname <- paste0("\\b",coef.report[j,1], "\\b")
      coef.report$Est[j] <- fixedef[grepl(fixed_names, pattern = parname),1]
      coef.report$SE[j] <- fixedef[grepl(fixed_names, pattern = parname),2]
    }
    
    # extract sigmas
    for (j in c(8,10,11,13,15,16,18)) {
      parname <- paste0("\\b",coef.report[j,1], "\\b")
      coef.report$Est[j] <- report[grepl(rep_names, pattern = parname),1]
      coef.report$SE[j] <- report[grepl(rep_names, pattern = parname),2]
    }
    print(paste0("model ", i))
    print(coef.report)
  }
}
