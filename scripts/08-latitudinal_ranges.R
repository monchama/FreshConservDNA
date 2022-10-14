library(readr)
library(visreg)
library(ggplot2)
library(rworldmap)
library(rgdal)
library(raster)
library(rgeos)
library(sp)
library(viridis)
library(dplyr)
library(ggthemes)
library(grid)
library(dismo)
library(plyr)
library(gridExtra)
library(colorspace)


#
raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


####bold data latitude range filling for all species
#estimate records pr latitude ban and compare to median range center

bold<-as.data.frame(read_tsv(file = file.choose())) #bold_records

#select Ca records
bold_ca<-bold[which(bold$country=="Canada"),]

bold_ca<-bold_ca[complete.cases(bold_ca[,48]),]


#create species list
ca_sp<-data.frame(table(bold_ca$species_name))

ca_sp$lat_span<-"NA"

ca_sp$cent_lat<-"NA"

ca_sp$phyllum<-"NA"

ca_sp$class<-"NA"


#calculate latitude span and centroid

for(i in 1:nrow(ca_sp)){
  
  ca_sp[i,3]<-max(bold_ca[which(bold_ca$species_name==ca_sp[i,1]),][,c(48)])-min(bold_ca[which(bold_ca$species_name==ca_sp[i,1]),][,c(48)])
  
  ca_sp[i,4]<- mean( as.matrix(unique(bold_ca[which(bold_ca$species_name==ca_sp[i,1]),][,c(48)])))
 
  ca_sp[i,5]<-bold_ca[which(bold_ca$species_name==ca_sp[i,1]),]$phylum_name[1]
  ca_sp[i,6]<-bold_ca[which(bold_ca$species_name==ca_sp[i,1]),]$class_name[1]
  
}

ca_sp$lat_span<- as.numeric(ca_sp$lat_span)

ca_sp$cent_lat<- as.numeric(ca_sp$cent_lat)


#nr of records relative to latitude span
ca_sp$dens<-ca_sp$Freq/(ca_sp$lat_span+0.1)

#plot
a1<-ggplot(ca_sp, aes(cent_lat, log(dens))) +
  geom_point(alpha=0.5) + 
  geom_vline(xintercept=57.5)+
  geom_hline(yintercept=1.75)+
  scale_color_viridis(discrete=TRUE) +
  xlab("Center of latitude range")+
  ylab("log(nr of bold records per latitude band)")+
  annotate("text", x = 71, y = -2.1, label = "Northern range",size = 4)+
  annotate("text", x = 45, y = -2.5, label = "Southern range",size = 4)+
  annotate("text", x = 70, y = 6, label = "High density",size = 4)+
  annotate("text", x = 70, y = -2.5, label = "Low density",size = 4)+
  ylim(-2.5,6)+
  xlim(40,75)+
  theme_classic()

a2<-ggplot(ca_sp, aes(cent_lat, log(dens),col=phyllum)) +
  geom_smooth(method=lm,se=F,lwd=2) + 
  scale_color_viridis(discrete=TRUE) +
  xlab("Center of latitude range")+
  ylab("log(nr of bold records per latitude band)")+
  ylim(-2.5,6)+
  xlim(40,75)+
  theme_classic()+
  theme(legend.position = c(0.8, 0.8))



grid.arrange(a1,a2,ncol=2)

##test for phyllum effect

m<-lm(log(dens)~cent_lat*phyllum,data=ca_sp)

plot(m)# quick inspection of the model residuals looks ok.
drop1(m,test="F")# no interaction effect


#looks like northern species has a lower density of records pr. latitude span


###############gbif work

###compare the distribution of bold and gbif records

#dowload a shapefile of Canada
worldMap <- getMap()

c_map <- worldMap[which(worldMap$NAME%in%c("Canada")),]

c_map<-spTransform(c_map, CRS("+proj=utm +zone=14 +datum=WGS84"))

plot(c_map)

#create a grid covering the spatial object
grid <- raster(extent(c_map), resolution = c(100000,100000), crs = proj4string(c_map))


# convert to SpatialPolygonsDataFrame
gridPolygon <- rasterToPolygons(grid)

#check if this looks ok
plot(c_map)
plot(gridPolygon, add = T)

#Clip the grid to fit the polygon boundary
intersectGridClipped <- raster::intersect(gridPolygon, c_map)

#check plot
plot(intersectGridClipped)

#gbif data
gbif<-read_tsv(file.choose()) #should be file containing all gbif obs


sp_gbif<-as.data.frame(gbif[,22:23])
sp_gbif<-sp_gbif[complete.cases(sp_gbif),]

coordinates(sp_gbif)<-cbind(sp_gbif$decimalLongitude,sp_gbif$decimalLatitude) 

sp_gbif<-sp_gbif[-which(sp_gbif$decimalLatitude==0),]
sp_gbif<-sp_gbif[-1229961,]

#extract bold Odonata records

sp_bold<- bold_ca[,48:49]
sp_bold<-sp_bold[complete.cases(sp_bold),]
coordinates(sp_bold)<-cbind(sp_bold$lon,sp_bold$lat) 


proj4string(sp_bold) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
sp_bold<-spTransform(sp_bold, CRS("+proj=utm +zone=14 +datum=WGS84"))


proj4string(sp_gbif) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
sp_gbif<-spTransform(sp_gbif, CRS("+proj=utm +zone=14 +datum=WGS84"))


intersectGridClipped$cellid<-1:nrow(intersectGridClipped)


#count the number of observation in each grid cell

bold_count <- over(sp_bold, intersectGridClipped)

tbold_count<-as.data.frame(table(bold_count$cellid))
colnames(tbold_count)<-c("cellid","bold_rec")


#########do spatial count for large dataset
split<-split(sp_gbif, (0:nrow(sp_gbif) %/% 100000)  )

totgbif<-as.data.frame(intersectGridClipped[,51])
totgbif$gbif_rec<-0


for(i in 1:length(split)){
gbif_count <- over(split[[i]], intersectGridClipped)
tgbif_count<-as.data.frame(table(gbif_count$cellid))
colnames(tgbif_count)<-c("cellid","gbif_rec")

totgbif<-rbind(totgbif,tgbif_count)
}

tgbif_count<- ddply(totgbif, ~cellid, summarise, gbif_rec = sum(gbif_rec))


######merge data counts

m1<-merge(intersectGridClipped, tgbif_count, by="cellid", all = TRUE)
m1<-merge(m1, tbold_count, by="cellid", all = TRUE)

m1<-as.data.frame(m1)
m1[is.na(m1$bold_rec),]$bold_rec<-0
m1[is.na(m1$gbif_rec),]$gbif_rec<-0

#plot relationship on log-log scale
ggplot(m1, aes(log(bold_rec+1), log(gbif_rec+1))) +
  geom_point() +
  geom_smooth(method = lm,se=F) + theme_classic()



#####################################
#########build a bivariate plot of GBIF and BOLD counts

#split each variable into 3 classes

m1$gbif_trans<-log(m1$gbif_rec+1)
m1$bold_trans<-log(m1$bold_rec+1)
m1$gbif_log<-log(m1$gbif_rec)
m1$bold_log<-log(m1$bold_rec)

h.v<-quantile(m1[which(m1$gbif_trans > 0),]$gbif_trans,c(0.5,0.75,1),na.rm=TRUE)
p.v<-quantile(m1[which(m1$bold_trans > 0),]$bold_trans,c(0.33,0.66,1),na.rm=TRUE)

m1<- m1 %>% mutate(y= ifelse(bold_trans<p.v[1],1,ifelse(bold_trans<p.v[2],2,3)) ,
                   x= ifelse(gbif_trans<h.v[1],1,ifelse(gbif_trans<h.v[2],2,3)))  


#overview plot of the splits
ggplot(data=m1,aes(x=gbif_trans,y=bold_trans,color=atan(y/x),alpha=x+y))+
  geom_point(size=2)+  
  guides(alpha="none",color="none")+
  geom_hline(yintercept=p.v,color="gray20",linetype=2)+
  geom_vline(xintercept=h.v,color="gray20",linetype=2)+
  scale_color_viridis(name="Color scale",option="magma")+theme_minimal()+
  labs(x="GBIF counts (log + 1 scale)",
       y="BOLD counts (log + 1 scale)")


#create final map
cagrid <- fortify(intersectGridClipped,region="cellid")

m1$id<-m1$cellid

cagrid  <-join(cagrid ,as.data.frame(m1), by="id")

cagrid$resid<-resid(lm(cagrid$gbif_trans~ cagrid$bold_trans))

###create dummy legend for plot
d<-expand.grid(x=1:3,y=1:3)
d<-merge(d,data.frame(x=1:3,xlabel=c(" ", " "," ")),by="x")
d<-merge(d,data.frame(y=1:3,ylabel=c(" ", " "," ")),by="y")

g.legend<-
  ggplot(d, aes(x,y,fill=atan(y/x),alpha=x+y,label=paste0(xlabel,"\n",ylabel)))+
  geom_tile()+
  geom_text(alpha=1)+
  scale_fill_viridis()+
  theme_void()+
  theme(legend.position="none",
        panel.background=element_blank(),
        plot.margin=margin(t=10,b=10,l=10))+
  labs(x="GBIF records (low to high)",y="BOLD records (low to high)")+
  theme(axis.title=element_text(color="black"))+
  # Draw some arrows:
  geom_segment(aes(x=1, xend = 3 , y=0, yend = 0), size=1.5,
               arrow = arrow(length = unit(0.6,"cm"))) +
  geom_segment(aes(x=0, xend = 0 , y=1, yend = 3), size=1.5,
               arrow = arrow(length = unit(0.6,"cm"))) 



#create bivariate map
bmap<-ggplot(data=cagrid , aes(y=lat, x=long, group=id, fill=atan(y/x),alpha=x+y)) +
  scale_fill_viridis_c() +
  labs(title="Bivariate map")+
  geom_polygon()+
  theme_map(base_size = 12) +
  theme(plot.title=element_text(size = 16, face="bold",margin=margin(b=10))) +
  theme(plot.subtitle=element_text(size = 14, margin=margin(b=-20))) +
  theme(plot.caption=element_text(size = 9, margin=margin(t=-15),hjust=0))  +  
  guides(alpha="none",fill="none")



######plot 3 panel

mapg<-ggplot(data=cagrid , aes(y=lat, x=long, group=id, fill=gbif_log)) +
  scale_fill_viridis_c(direction = -1) +
  labs(title="GBIF")+
  geom_polygon()+
  theme_map(base_size = 12)+
  theme(legend.position = c(0.8, 0.5),legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'))

mapb<-ggplot(data=cagrid , aes(y=lat, x=long, group=id, fill=bold_log)) +
  scale_fill_viridis_c(direction = -1) +
  labs(title="BOLD")+
  geom_polygon()+
  theme_map(base_size = 12) +
  theme(legend.position = c(0.8, 0.5),legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.2, 'cm'))


grid.arrange(mapg,mapb,bmap)

vp<-viewport(width=0.5,height=0.15,x=0.68,y=0.24)

print(g.legend+labs(title=""),vp=vp)
