#Are there any species that possibly went extinct?

data<-readRDS("data/compiled_data.RDS")

str(data)

last_observed <- data %>%
  filter(count>0)

last_observed <- last_observed %>%
  group_by(para.spc) %>%
  summarize(max_year = max(year))

last_observed_1980 <- last_observed %>%
  filter(max_year<1981)
last_observed_1980<-as.data.frame(last_observed_1980)
last_observed_1980<-last_observed_1980[order(last_observed_1980$max_year),]

last_observed_1970 <- last_observed %>%
  filter(max_year<1971)

last_observed_1950<- last_observed %>%
  filter(max_year<1951)

thelist<-last_observed_1980$para.spc

filtered_data <- data %>%
  filter(para.spc %in% thelist)

stuffs <- filtered_data %>%
  mutate(occurrence = (count>0))
stuffs$occur <- as.numeric(stuffs$occurrence)


# To make a plot of occurrence~year

rep_data <- filtered_data %>%
  group_by(para.spc) %>%
  summarize(replicates = n())

sp_1_data <- stuffs %>%
  filter(para.spc == thelist[1])

sp_1_tally <- filtered_data %>%
  filter(para.spc == thelist[1])

sp_1_tally$before_after <- if_else(sp_1_tally$year > last_observed_1980$max_year[1], "after", "before")

tally <- sp_1_tally %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_2_data <- stuffs %>%
  filter(para.spc == thelist[2])

sp_2_data$before_after <- if_else(sp_2_data$year > last_observed_1980$max_year[2], "after", "before")

tally <- sp_2_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_3_data <- stuffs %>%
  filter(para.spc == thelist[3])

sp_3_data$before_after <- if_else(sp_3_data$year > last_observed_1980$max_year[3], "after", "before")

tally <- sp_3_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_4_data <- stuffs %>%
  filter(para.spc == thelist[4])

sp_4_data$before_after <- if_else(sp_4_data$year > last_observed_1980$max_year[4], "after", "before")

tally <- sp_4_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_5_data <- stuffs %>%
  filter(para.spc == thelist[5])

sp_5_data$before_after <- if_else(sp_5_data$year > last_observed_1980$max_year[5], "after", "before")

tally <- sp_5_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_6_data <- stuffs %>%
  filter(para.spc == thelist[6])

sp_6_data$before_after <- if_else(sp_6_data$year > last_observed_1980$max_year[6], "after", "before")

tally <- sp_6_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_7_data <- stuffs %>%
  filter(para.spc == thelist[7])

sp_7_data$before_after <- if_else(sp_7_data$year > last_observed_1980$max_year[7], "after", "before")

tally <- sp_7_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_8_data <- stuffs %>%
  filter(para.spc == thelist[8])

sp_8_data$before_after <- if_else(sp_8_data$year > last_observed_1980$max_year[8], "after", "before")

tally <- sp_8_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_9_data <- stuffs %>%
  filter(para.spc == thelist[9])

sp_9_data$before_after <- if_else(sp_9_data$year > last_observed_1980$max_year[9], "after", "before")

tally <- sp_9_data %>%
  group_by(before_after) %>%
  summarize(count = n())

sp_10_data <- stuffs %>%
  filter(para.spc == thelist[10])

sp_10_data$before_after <- if_else(sp_10_data$year > last_observed_1980$max_year[10], "after", "before")

tally <- sp_10_data %>%
  group_by(before_after) %>%
  summarize(count = n())


library(wesanderson)
pal<-wes_palette(name="Zissou1",10,type="continuous")
str(pal)

thelist[1]
sp_1_plot<-ggplot(data = sp_1_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("COP-HP-LEP")+
  geom_vline(xintercept = last_observed_1980$max_year[1],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[1], color="#EA5C00", x = (last_observed_1980$max_year[1]+2), y = 0.5, hjust=0)+
  #annotate("text",label=rep_data$replicates[1],x=2010,y=1,hjust=0,size=10)+
  #ggtitle("Lepeophtheirus sp., an adult copepod in H. pretiosus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#E4B80E",0.3)),
        panel.background = element_rect(fill = alpha("#E4B80E",0.02),colour = alpha("#E4B80E",0.02)))

thelist[2]
sp_2_plot<-ggplot(data = sp_2_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("CES-HP-SP1")+
  geom_vline(xintercept = last_observed_1980$max_year[2],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[2], color="#EA5C00", x = (last_observed_1980$max_year[2]+2), y = 0.5, hjust=0)+
  #ggtitle("cestode sp. 1, a larval cestode in H. pretiosus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[3]
sp_3_plot<-ggplot(data = sp_3_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("TRE-SC-OPE")+
  geom_vline(xintercept = last_observed_1980$max_year[3],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[3], color="#EA5C00", x = (last_observed_1980$max_year[3]+2), y = 0.5, hjust=0)+
  #ggtitle("Opechona sp., an adult trematode in S. caurinus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[4]
sp_4_plot<-ggplot(data = sp_4_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("TRE-HC-GON")+
  geom_vline(xintercept = last_observed_1980$max_year[4],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[4], color="#EA5C00", x = (last_observed_1980$max_year[4]+2), y = 0.5, hjust=0)+
  #ggtitle("Gonocerca sp., an adult trematode in H. colliei")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[5]
sp_5_plot<-ggplot(data = sp_5_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("TRE-HP-MET")+
  geom_vline(xintercept = last_observed_1980$max_year[5],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[5], color="#EA5C00", x = (last_observed_1980$max_year[5]+2), y = 0.5, hjust=0)+
  #ggtitle("metacercaria sp., a larval trematode in H. pretiosus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[6]
sp_6_plot<-ggplot(data = sp_6_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("CES-MP-NYB")+
  geom_vline(xintercept = last_observed_1980$max_year[6],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[6], color="#EA5C00", x = (last_observed_1980$max_year[6]+2), y = 0.5, hjust=0)+
  #ggtitle("Nybelinia surmenicola, a larval cestode in M. productus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.6)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.01),colour = alpha("#3B9AB2",0.01)))

thelist[7]
sp_7_plot<-ggplot(data = sp_7_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("TRE-SC-SP2")+
  geom_vline(xintercept = last_observed_1980$max_year[7],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[7], color="#EA5C00", x = (last_observed_1980$max_year[7]+2), y = 0.5, hjust=0)+
  #ggtitle("trematode sp. 2, an adult trematode in S. caurinus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[8]
sp_8_plot<-ggplot(data = sp_8_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(labels=NULL, limits = c(1880,2020))+
  xlab("")+
  ylab("TRE-PV-SP1")+
  geom_vline(xintercept = last_observed_1980$max_year[8],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[8], color="#EA5C00", x = (last_observed_1980$max_year[8]+2), y = 0.5, hjust=0)+
  #ggtitle("trematode sp. 1, an adult trematode in P. vetulus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[9]
sp_9_plot<-ggplot(data = sp_9_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  xlim(1880,2020)+
  xlab("")+
  ylab("TRE-SC-MET")+
  geom_vline(xintercept = last_observed_1980$max_year[9],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[9], color="#EA5C00", x = (last_observed_1980$max_year[9]+2), y = 0.5, hjust=0)+
  #ggtitle("metacercaria sp., a larval trematode in S. caurinus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

thelist[10]
sp_10_plot<-ggplot(data = sp_10_data, aes(x = year, y = jitter(occur,0.2)))+
  geom_point()+
  scale_y_continuous(breaks=c(0,1))+
  xlim(1880,2020)+
  xlab("")+
  ylab("TRE-SC-SP1")+
  geom_vline(xintercept = last_observed_1980$max_year[10],lty=2, color="#EA5C00")+
  annotate("text",label=last_observed_1980$max_year[10], color="#EA5C00", x = (last_observed_1980$max_year[10]+2), y = 0.5, hjust=0)+
  #ggtitle("trematode sp. 1, an adult trematode in S. caurinus")+
  theme_classic()+
  theme(plot.background = element_rect(fill = alpha("#3B9AB2",0.3)),
        panel.background = element_rect(fill = alpha("#3B9AB2",0.005),colour = alpha("#3B9AB2",0.005)))

library(cowplot)
extinct_panels <- ggdraw(plot=NULL,xlim=c(0,20),ylim=c(0,50))+
  draw_plot(sp_1_plot,x=1.25,y=40.5,width=9,height=9)+
  draw_plot(sp_2_plot,x=10.5,y=40.5,width=9,height=9)+
  draw_plot(sp_3_plot,x=1.25,y=31,width=9,height=9)+
  draw_plot(sp_4_plot,x=10.5,y=31,width=9,height=9)+
  draw_plot(sp_5_plot,x=1.25,y=21.5,width=9,height=9)+
  draw_plot(sp_6_plot,x=10.5,y=21.5,width=9,height=9)+
  draw_plot(sp_7_plot,x=1.25,y=12,width=9,height=9)+
  draw_plot(sp_8_plot,x=10.5,y=12,width=9,height=9)+
  draw_plot(sp_9_plot,x=1.25,y=2.5,width=9,height=9)+
  draw_plot(sp_10_plot,x=10.5,y=2.5,width=9,height=9)+
  draw_label("year",x=10,y=1.5,size=30)+
  draw_label("occurrence",x=0.5,y=25,size=30,angle=90)
extinct_panels

