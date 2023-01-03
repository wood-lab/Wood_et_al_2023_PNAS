### Figure 2
### Created by Chelsea Wood
### chelwood@uw.edu

library(tidyverse)
library(cowplot)

# Load the parameter estimates in the file phase_1_results.Rdata, which creates a list object.

indiv_psite_data<-read.csv("data/for_plotting/spp_level_rand_effects.csv",sep=",",header=T)
psite_type_data<-read.csv("data/corresponding_psite_types_tim_to_universal_code.csv",sep=",", header=T)
psite_host_data<-read.csv("data/psite_lc_complexity.csv",sep=",", header=T)
final_data<-merge(psite_type_data,indiv_psite_data,by.x="tim_code",by.y="X",all=TRUE)
final_data<-merge(final_data,psite_host_data,by.x="tim_code",by.y="tim_code",all=TRUE)
length(final_data$tim_code)
final_data<-final_data[,-c(7,8)]
names(final_data)<-c("tim_code","psite_type","psite_code","label","Estimate","Std..Error","n_hosts","terminal")


# Do a quick tally

tallies <- final_data %>%
  group_by(psite_type) %>%
  count()
tallies


# Put everything in order so you don't have to do it manually later on

final_data_ordered <- final_data %>%
  arrange(psite_code) %>%
  arrange(factor(n_hosts, levels = c("1","2","3")))
final_data_ordered$order<-rev(1:85)

# Get your colors

library(wesanderson)
pal<-wes_palette(name="Zissou1",10,type="continuous")
cols <- wes_palette(10, name = "Zissou1",type="continuous")[c(1,4,7)]


# Organized by number of hosts

indiv_psites_plot<-ggplot(final_data_ordered,aes(x=Estimate,y=psite_code,label=n_hosts))+
  geom_rect(xmin=-2.09,xmax=2.09,ymin=65.5,ymax=85.6,fill="#E4B80E",alpha=0.05)+
  geom_rect(xmin=-2.09,xmax=2.09,ymin=44.5,ymax=65.5,fill="#9EBE91",alpha=0.05)+
  geom_rect(xmin=-2.09,xmax=2.09,ymin=0.4,ymax=44.5,fill="#3B9AB2",alpha=0.05)+
  geom_point()+
  geom_errorbar(aes(xmin=Estimate-Std..Error,xmax=Estimate+Std..Error))+
  geom_vline(xintercept = 0,lty=1)+
  geom_hline(yintercept = 65.5,lty=3,lwd=0.25)+
  geom_hline(yintercept = 44.5,lty=3,lwd=0.25)+
  xlab("random effect of year")+
  ylab("parasite code")+
  xlim(-1.9,1.9)+
  scale_y_discrete(limits=rev(final_data_ordered$psite_code))+
  #geom_text(x=-2,angle=0,hjust=0,vjust=0.5)+
  #annotate("text",label=c("Copepoda","Hirudinea","Monogenea","Trematoda","Cestoda","Nematoda","Acanthocephala"),
  #         x=,y=(final_data_ordered$order+1.35),hjust=0.5,vjust=0.5,size=3)+
  annotate("text",label=c("1 host","2 hosts","3+ hosts"),x=-1.8,y=c(76,54,23),size=5,hjust = 0.5)+
  theme_classic()+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),panel.grid.major.y = element_line(color="darkgray"),
        axis.text.y = element_text(size=7), axis.title = element_text(size = 18))+
  coord_cartesian(clip="off")

indiv_psites_plot

