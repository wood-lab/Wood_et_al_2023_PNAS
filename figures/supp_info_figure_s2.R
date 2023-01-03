### Supplementary Information Figure S2
### Created by Tim Essington
### essing@uw.edu

library(TMB)
library(ggplot2)
library(tidyr)
library(viridis)
# load phase 2 results data file
load("src/analysis/phase_2_results.Rdata")

year.cont.model <- models$cont$contyear 
re <- summary(year.cont.model$rep, "random")
g1 <- re[grep(rownames(re), pattern = "\\bg1\\b"),]
g2 <- re[grep(rownames(re), pattern = "\\bg3\\b"),]
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
  xlab(bquote(gamma[k*","*cont_ord1])) +
  ylab(bquote(gamma[k*","*year])) +
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())

