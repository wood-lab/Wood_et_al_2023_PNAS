### Figure 1
### Created by Chelsea Wood
### chelwood@uw.edu

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(rgdal)
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(ggmap)
library(ggsn)
library(dplyr)
library(ggmap)
library(cowplot)
library(sp)
library(rgeos)
library(raster)
library(maptools)
library(rgdal)
library(wesanderson)
library(dplyr)
library(cowplot)


#Read in data

lots<-read.csv("data/for_plotting/all_fish_spp.csv")
str(lots)

lots$Species.name<-gsub("Clupea palasii pallasii","Clupea pallasii",lots$Family.Species)
lots$Species.name<-gsub("Clupea palasii ","Clupea pallasii",lots$Species.name)
lots$Species.name<-gsub("Theragra chalcogramma", "Gadus chalcogrammus",lots$Species.name)

#Find the right colors

wes_palette("Zissou1")
pal<-wes_palette(name = "Zissou1", 9, type = "continuous")

#Map

bounds<-c(left=-123.3, bottom=47, right=-122.1, top=49.0)
map<-get_stamenmap(bounds, zoom=8, maptype = "terrain-background") %>% ggmap()+
  geom_point(data = lots, aes(x=jitter(long,factor=1000),y=jitter(lat,factor=100),fill=Year.Collected),size=2.5,shape=21,
             alpha=0.7)+
  scale_fill_gradient("Year",low="black",high="white")+
  xlab("Longitude (°W)")+
  ylab("Latitude (°N)")+
  theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text=element_text(size=14),axis.title=element_text(size=14))
  
map


#Reference map

bounds<-c(left=-125, bottom=25, right=-66.6, top=50)
ref_map<-get_stamenmap(bounds, zoom=5, maptype = "terrain-background") %>% ggmap()+
  geom_rect(xmin=-123.3,ymin=47,xmax=-122.1,ymax=49,color="black",fill=NA)+
  scale_fill_gradient("Year",low="black",high="white")+
  xlab("")+
  ylab("")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text=element_blank(),axis.title=element_blank(),axis.ticks = element_blank())
ref_map


#Plot

library(viridis)
viridis(8)


plot<-ggplot(data=lots,aes(x=Year.Collected,y=lat))+
  geom_point(data=lots,aes(x=jitter(Year.Collected,factor=30),y=jitter(lat,factor=100),fill=Species.name),shape=21,size=4)+
  scale_fill_manual(limits=c("Clupea pallasii","Embiotoca lateralis","Gadus chalcogrammus",
                             "Hydrolagus colliei","Hypomesus pretiosus",
                             "Merluccius productus","Parophrys vetulus","Sebastes caurinus")
                    ,labels=c("C. pallasii","E. lateralis","G. chalcogrammus",
                              "H. colliei","H. pretiosus",
                              "M. productus","P. vetulus","S. caurinus")
                    ,values=c("#440154FF", "#46337EFF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF"),
                    name="Host spp.")+
  xlab("")+
  ylab("Latitude (°N)")+
  ylim(47,49)+
  xlim(1879,2020)+
  theme_bw()+
  theme(legend.position = "right",legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(face = "italic"))+
  guides(fill = guide_legend(override.aes = list(size=7)))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
plot


# Now make the environmental variables plot

# Read in the temp data

RR_temp<-read.csv("data/for_plotting/RR_temp.csv")
# Need to scale all of the variables between zero and one

RR_temp$sTemp <- 1/(max(RR_temp$temp_na_rm, na.rm=TRUE)-min(RR_temp$temp_na_rm, na.rm=TRUE))*RR_temp$temp_na_rm-
  min(RR_temp$temp_na_rm, na.rm=TRUE)/(max(RR_temp$temp_na_rm, na.rm=TRUE)-min(RR_temp$temp_na_rm, na.rm=TRUE))


#Pollutant data

compiled_data <-readRDS("data/compiled_data.RDS")

pollutants <- compiled_data %>%
  group_by(year) %>%
  summarise(ord1_annual = mean(ord1),
            ord2_annual = mean(ord2))

pollutants$sOrd1 <- (pollutants$ord1_annual+abs(min(pollutants$ord1_annual, na.rm=TRUE)))/max(pollutants$ord1_annual+abs(min(pollutants$ord1_annual, na.rm=TRUE)),na.rm=T)
pollutants$sOrd2 <- (pollutants$ord2_annual+abs(min(pollutants$ord2_annual, na.rm=TRUE)))/max(pollutants$ord2_annual+abs(min(pollutants$ord2_annual, na.rm=TRUE)),na.rm=T)

#Stole a bunch of code from Tim Essington (whose Figure 1 inspired the figure I want to make)

require(readxl)
require(dplyr)

# function to interpolate in between data years for 'nsteps'.  Makes plot look smoother
interp.xy<-function(x,y,inc=0.1){ 
  xout <- seq(min(x), max(x), by = inc)
  interpolated <- approx(x, y, xout)
  return(list(x=interpolated$x,y=interpolated$y))
}

#base function for plotting.  It interpolates given time series of x and y, makes 
# a polygon on an existing plot.  The parameter 'mult' essentially scales width of the polygons,
# and the parameter base base refers to the position of the polygon on the y axis.

plot.cartoon <- function(base, data, color.list, y.list, mult=0.3) {
  
  x <- as.array(as.matrix(data[,1]))
  y <- as.array(as.matrix(data[,2]))
  new.xy<-interp.xy(x=x,y=y) 
  tmp.data<-cbind(new.xy$x,new.xy$y)
  n.points<-length(new.xy$x)
  # makes a tiny polygon for each increment, setting the polygon width
  # and color equal to the vale of the metric
  for (i in 2:n.points){
    poly.col <- color.list[round(new.xy$y[i], digits = 2)  == y.list] # look up color to use based on value of the y variable
    xplot<-c(rep(tmp.data[i-1,1],2),rep(tmp.data[i,1],2)) # create x values for polygon
    yplot<-c(-mult*tmp.data[i,2],mult*tmp.data[i,2], mult*tmp.data[i,2],-mult*tmp.data[i,2])+rep(base,4) # create y values of polygon, where mult just scales the width of the polygon
    polygon(x = xplot,
            y = yplot,
            border = poly.col,
            col = poly.col)
  }
}

# set up parameters for plotting
y.list <- round(seq(0, 1, by = 0.01), digits = 2)
#col.pallete<-colorRamp(c("blue","yellow","green","orange","red"),space="rgb",interpolate="linear",bias=1) # this is a color version that we removed
color.palette<-colorRampPalette(c("grey90","black"))
color.list <- color.palette(length(y.list))

# compile the data to use in plot
par(mai=c(1,2.5,0,0.5))
# create empty plot to insert polygons
plot(x=0, y=0, type="n", xaxs="i", yaxs="i", las=1, 
     xlim=c(1880,2020), ylim=c(0,10),xlab="Year",axes=FALSE,ylab="",cex.lab=2.5)
axis(1,at=seq(1880,2020,5),cex.axis=2)
par(xpd=TRUE)

base.list<-seq(1,16,1)

text(x=1880,y=base.list[9],labels="Annual temperature",pos=2,cex=1.5)
text(x=1880,y=base.list[8],labels="Pollutant PC1",pos=2,cex=1.5)
text(x=1880,y=base.list[7],labels="Pollutant PC2",pos=2,cex=1.5)
text(x=1880,y=base.list[5],labels=substitute(paste(italic("G. chalcogrammus"))),pos=2,cex=1.5)
text(x=1880,y=base.list[1],labels=substitute(paste(italic("P. vetulus"))),pos=2,cex=1.5)
text(x=1880,y=base.list[2],labels=substitute(paste(italic("M. productus"))),pos=2,cex=1.5)
text(x=1880,y=base.list[3],labels=substitute(paste(italic("H. pretiosus"))),pos=2,cex=1.5)
text(x=1880,y=base.list[4],labels=substitute(paste(italic("H. colliei"))),pos=2,cex=1.5)
text(x=1880,y=base.list[6],labels=substitute(paste(italic("C. pallasii"))),pos=2,cex=1.5)


data.2.use <- c("sOrd1","sOrd2")
base.list<-c(8,7) # the list of polygon centers on the y axis for each metric
for (j in 1:2){
  base<-base.list[j]
  data <- pollutants %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sTemp")
base.list<-c(9) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- RR_temp %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(YEAR, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

# Bring in fish density data

herring<-read.csv("data/processed/for_plotting/herring_density.csv",sep=",",header=TRUE)
herring$sHerring <- herring$fish_density/max(herring$fish_density, na.rm=TRUE)
smelt<-read.csv("data/processed/for_plotting/smelt_density.csv",sep=",",header=TRUE)
smelt$sSmelt <- smelt$fish_density/max(smelt$fish_density, na.rm=TRUE)
hake<-read.csv("data/processed/for_plotting/hake_density.csv",sep=",",header=TRUE)
hake$sHake <- hake$fish_density/max(hake$fish_density, na.rm=TRUE)
pollock<-read.csv("data/processed/for_plotting/pollock_density.csv",sep=",",header=TRUE)
pollock$sPollock <- pollock$fish_density/max(pollock$fish_density, na.rm=TRUE)
ratfish<-read.csv("data/processed/for_plotting/ratfish_density.csv",sep=",",header=TRUE)
ratfish$sRatfish <- ratfish$fish_density/max(ratfish$fish_density, na.rm=TRUE)
sole<-read.csv("data/processed/for_plotting/sole_density.csv",sep=",",header=TRUE)
sole$sSole <- sole$fish_density/max(sole$fish_density, na.rm=TRUE)

# Add the species in this order:
# Clupea - herring
# Hydrolagus - ratfish
# Hypomesus - smelt
# Merluccius - hake
# Parophrys - sole
# Theragra - pollock

data.2.use <- c("sHerring")
base.list<-c(6) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- herring %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sRatfish")
base.list<-c(4) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- ratfish %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sSmelt")
base.list<-c(3) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- smelt %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sHake")
base.list<-c(2) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- hake %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sSole")
base.list<-c(1) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- sole %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

data.2.use <- c("sPollock")
base.list<-c(5) # the list of polygon centers on the y axis for each metric
for (j in 1:1){
  base<-base.list[j]
  data <- pollock %>%
    filter (data.2.use[j]!="NA") %>%
    dplyr::select(year, data.2.use[j])
  plot.cartoon(base, data, color.list, y.list = y.list)
  
}

par(xpd=TRUE)


# Now put everything together into a single, four-panel figure

final_figure <- ggdraw(plot=NULL,xlim=c(0,20),ylim=c(0,20))+
  draw_image("figures/env_plot_full_v7.jpeg",x=6.75,y=-1.85,width=13,height=13)+
  draw_plot(map,x=-1.35,y=-0.18,width=10,height=20)+
  draw_plot(ref_map,x=4.5,y=15,width=3,height=3)+
  draw_plot(plot,x=8.1,y=9.75,width=12,height=10)+
  draw_label("(a)",x=0.4,y=19.5,size=30)+
  draw_label("(b)",x=7.8,y=19.5,size=30)+
  draw_label("(c)",x=7.8,y=10,size=30)
final_figure

dev.off()
