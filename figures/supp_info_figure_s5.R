### Supplementary Information Figure S5
### Created by Tim Essington
### essing@uw.edu

library(dplyr)

source("src/load_data.R")

spc <- levels(thedata$para.spc)
nspc <- length(spc)
ngrp <- ncol(A)

# create a vector of parasite species ID for each observation
Xspc <- (as.numeric(thedata$para.spc) -1)
spc.groups <- (levels(thedata$par.type))

# get index of taxonomic group for each species, and create species specific LH type
# and number of hosts in life cycle
NH.spc.tmp <- rep(NA, nspc)
for (i in 1:length(spc)) {
  NH.spc.tmp[i] <- dplyr::filter(thedata, para.spc == spc[i])$nhosts[1]
}

# subtract 1 for C++ indexing

NHspc <- as.integer(NH.spc.tmp -1)
load("src/analysis/phase_1_results.Rdata")

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
                            
