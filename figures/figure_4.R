### Figure 4
### Created by Tim Essington and Chelsea Wood
### essing@uw.edu and chelwood@uw.edu

library(tidyverse)
library(ggplot2)
library(cowplot)
library(viridis)
library(gridExtra)

# Start by making a figure that summarizes the slopes of the various putative drivers

# Load the parameter estimates in the file phase_2_results.Rdata, which creates a list object.
#load("~/Dropbox/Vault/University of Washington/Projects/Historical parasites/Innovation Award/PoP_Allspp/src/analysis/phase_2_results.Rdata")

# now load phase 2 results
load("src/analysis/phase_2_results.Rdata")
driver_slopes <- models$fit_summary

# That will give you a dataframe that has the estimated coefficient for each predictor variable, and its SE.  
# Note that these have unique estimates for each taxonomic group, so you'll probably want to pull out the mean:

drivers_plot<-ggplot(driver_slopes,aes(x=predictor,y=slope))+
  geom_point(size=4)+
  geom_errorbar(aes(ymin=slope-se,ymax=slope+se,width=0.3))+
  geom_hline(yintercept = 0,lty=1)+
  xlab("potential driver")+
  ylab("slope")+
  scale_x_discrete(limits=c("year","temperature","contaminant PC1","contaminant PC2","density"),labels=c("year"="year","temp"="temp",
                                                                                      "contaminant PC1" = "PC1",
                                                                                      "contaminant PC2" = "PC2",
                                                                                      "density" = "host density"))+
  theme_classic()+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid = element_blank())+
  coord_cartesian(clip="off")

drivers_plot


# Then we'll make a plot that shows the bivariate relationship between the most important driver (temperature) and parasite
# burden.

load(file = "src/analysis/temp_predictions.Rdata")

# This will also include the reduced data set that only has complex life cycles and years with temperature data (red.data)
# If you try to plot it, it's a mess because the data are so highly dispersed

library(wesanderson)
pal<-wes_palette(name="Zissou1",8,type="continuous")
pal<-rev(pal)


# Try binning to reduce noise. Start with 0.1-degree increments of temperature, ignoring parasite type.

red.data$temp_round <- round(red.data$temp_na_rm, digits = 1)
red.data$year_round <- round(red.data$year, digits = -1)
binned <- red.data %>%
  group_by(temp_round) %>%
  summarize(bin_count = mean(count), bin_se = sd(count)/sqrt(n()), bin_n = n())

binned <- binned %>%
  filter(bin_n>25)


binned_year <- red.data %>%
  group_by(year_round) %>%
  summarize(bin_count = mean(count), bin_se = sd(count)/sqrt(n()), bin_n = n())

binned_year <- binned_year %>%
  filter(bin_n>25)

# now make the final plot
axis.scale <- .2
temp_plot<-ggplot(binned, aes(x = jitter(temp_round, factor=0), y = jitter(bin_count, factor=0))) +
  scale_y_continuous(limits = c(0, 8), expand = expansion(mult = c(0, 0.0)), oob = scales::squish,
                     sec.axis = sec_axis(~.*axis.scale, name="marginal effect", breaks = c(0,0.5, 1, 1.5, 2))) +
  geom_line(data = predict_df_temp, aes(x = temp, y = yhat / axis.scale), inherit.aes = FALSE, lwd = 1.5)+ylab("parasite count")+
  geom_ribbon(
    data = predict_df_temp,
    aes(
      x = temp,
      ymin = lb/axis.scale,
      ymax = ub/axis.scale
    ),
    inherit.aes = FALSE,
    fill = "darkgray",
    alpha = 0.4,
    show.legend = F
  )+
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = bin_count-bin_se, ymax = bin_count+bin_se))+
  xlab("SST (°C)")+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())
temp_plot

axis.scale.year <- 0.25
year_plot<-ggplot(binned_year, aes(x = jitter(year_round, factor=0), y = jitter(bin_count, factor=0))) +
  scale_y_continuous(limits = c(0, 4), expand = expansion(mult = c(0, 0.0)), oob = scales::squish,
                     sec.axis = sec_axis(~.*axis.scale.year, name="marginal effect", breaks = c(0,0.5, 1,1.5, 2))) +
  geom_line(data = predict_df_year, aes(x = year, y = yhat / axis.scale.year), inherit.aes = FALSE, lwd = 1.5)+ylab("parasite count")+
  geom_ribbon(
    data = predict_df_year,
    aes(
      x = year,
      ymin = lb/axis.scale.year,
      ymax = ub/axis.scale.year
    ),
    inherit.aes = FALSE,
    fill = "darkgray",
    alpha = 0.4,
    show.legend = F
  )+
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = bin_count-bin_se, ymax = bin_count+bin_se))+
  xlab("year")+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())
year_plot


# make 2-d plot of species slopes
temp.year.model<- models$temp$tempyear
re <- summary(temp.year.model$rep, "random")
g1 <- re[grep(rownames(re), pattern = "\\bg1\\b"),]
g2 <- re[grep(rownames(re), pattern = "\\bg2\\b"),]
## Prepare plot for analysis
gs <- tibble(g.t.hat = g1[,1],
             g.t.se = g1[,2],
             g.y.hat = g2[,1],
             g.y.se = g2[,2]
)
# get gbars
fe <- summary(temp.year.model$rep, "fixed")
gbar1 <- fe[grep(rownames(fe), pattern = "gbar1"),]
gbar2 <- fe[grep(rownames(fe), pattern = "gbar2"),]


gbars <- c(gbar1, gbar2)
names(gbars) <- c("gbar1", "gbar1se", "gbar2", "gbar2se")
gbars <- tibble(gbar1hat = gbars[1],
                gbar1se = gbars[2],
                gbar2hat = gbars[3],
                gbar2se = gbars[4])


col <- viridis(n = 10, option = "mako", alpha = 0.3) [3]
colse <- viridis(n = 10, option = "mako", alpha = 0.3) [7]
colgbar <- viridis(n = 10, option = "rocket") [3]
colgse <- viridis(n = 10, option = "rocket") [5]


species_plot <- ggplot(gs, aes(x = g.t.hat, y = g.y.hat)) +
  geom_errorbar(data = gs,aes(ymin = g.y.hat - g.y.se,
                    ymax = g.y.hat + g.y.se),
                size = 0.35,
                color = colse) + 
  geom_errorbar(data = gs,
                mapping = aes(xmin = g.t.hat - g.t.se,
                    xmax = g.t.hat + g.t.se),
                size = 0.35,
                color = colse) + 
  geom_point(colour = col, size = 1.5) + 
   ylim(-1.7, 1.25) + 
    xlim(-.5, 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_errorbar(data = gbars, mapping = aes(y = gbar2hat, xmin = gbar1hat - gbar1se,
                                  xmax = gbar1hat + gbar1se),
                size = 1,
                color = colgse, 
                inherit.aes = FALSE,
                width = 0) + 
  geom_errorbar(data = gbars, mapping = aes(x = gbar1hat, ymin = gbar2hat - gbar2se,
                                            ymax = gbar2hat + gbar2se),
                size = 1,
                color = colgse, 
                inherit.aes = FALSE,
                width = 0) + 
  geom_point(data = gbars, aes(x = gbar1hat, y = gbar2hat),
             colour = colgbar,
             size = 5.5) + 
  xlab(bquote(gamma[k*","*temp])) +
  ylab(bquote(gamma[k*","*year])) +
  theme_classic()+
  theme(legend.position = "top", legend.spacing.x = unit(0.001, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12), axis.title = element_text(size = 18),
        panel.grid.major.x = element_blank())
species_plot

# Now put it all together
png(file = "figures/phase_2_results.png",
    width = 9,
    height = 8,
    units = "in",
    res = 300
)

drivers_panels <- ggdraw(plot=NULL,xlim=c(0,30),ylim=c(0,10))+
  #draw_image("figures/fitted_time_effects.jpeg",x=-1.2,y=9.25,width=11.5,height=10,scale=1.05)+
  draw_plot(species_plot,x=0,y=0,width=8,height=10)+
  draw_plot(temp_plot,x=8,y=0,width=11,height=10)+
  draw_plot(year_plot,x=19,y=0,width=11,height=10)+
  draw_label("(a)",x=0.6,y=9.5,size=24)+
  draw_label("(b)",x=8.2,y=9.5,size=24) +
  draw_label("(c)",x=19.2,y=9.5,size=24)
  #draw_label("(d)",x=10.4,y=9.5,size=24)
drivers_panels
dev.off()
