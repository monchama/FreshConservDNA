############## State of biodiversity genomic records in Canada
######### Authors: Olivier Morissette, UQAC
### Last updated: July 25th 2022

rm(list=ls())
library(dplyr)
library(vegan)
library(sf)
library(raster)
library(terra)
library(fossil)
library(viridis)
library(ggplot2)
library(ggspatial)
library(gridExtra)
library(multcomp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Loading the raster layers ######
layers_landcover <- rast(list.files("../02_donnees/landcover/", pattern = "\\.tif$", full.names = TRUE))
layers_footprint <- rast(list.files("../02_donnees/human_footprint/", pattern = "\\.tif$", full.names = TRUE))
layers_can_threat <- rast(list.files("../02_donnees/canada_threat/", pattern = "\\.tif$", full.names = TRUE))

# Projecting the layer (this can take time) #
layers_p_landcover<-project(layers_landcover,"+proj=longlat +datum=WGS84", method="ngb")
layers_p_footprint<-project(layers_footprint,"+proj=longlat +datum=WGS84", method="ngb")
layers_p_threat<-project(layers_can_threat,"+proj=longlat +datum=WGS84", method="ngb")

######### Loading  and filtering the bold records #############
bold_records<-read.csv("../02_donnees/2022-06-18_bold_records.tsv", sep="\t", dec = ".", stringsAsFactors = T)
bold_records$coord<-is.na(bold_records$lon)
bold_records_occ<-filter(bold_records, coord == "FALSE")%>% droplevels() #to remove records without coordinates


###### Extracting land use and human footprint information ############
bold_records_project <- st_as_sf(bold_records_occ, coords = c("lon","lat")) #transforming the records to sf objects
bold_records_project <- st_set_crs(bold_records_project, 4326) #projecting to WGS-84

c<-terra::extract(layers_p_landcover, vect(bold_records_project), method="simple")
d<-terra::extract(layers_p_footprint, vect(bold_records_project), method="simple")
e<-terra::extract(layers_p_threat, vect(bold_records_project), method="simple")


bold_records_occ$land_use<-c[,2]
bold_records_occ$human_footprint<-d[,2]
bold_records_occ$canada_threat<-e[,2]


###### Writing the new dataset ########
write.csv2(bold_records_occ,"../02_donnees/2022-06-18_bold_records_environ_threat.csv")
rm(list=ls()) #cleaning the workspace


######## Boxplots for figure 2########

library(wesanderson) #for the color palette
data<-read.csv2("../02_donnees/2022-07-18_bold_records_environ_threat.csv", stringsAsFactors = T)
data$phylum_name<-factor(data$phylum_name, levels = c("Arthropoda","Magnoliophyta","Chordata","Mollusca","Pteridophyta",
                                                      "Lycopodiophyta","Pinophyta","Charophyta")) #ordering the factors to fit the symbology

###Boxplot of human footprint by phylum for world dataset###

pal8<-c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')

png("../04_figures/20220704_boxplot_footprint_global.png", width = 4,height = 7, units = "in", res=300)
par(mar=c(4,8.4,0.5,0.32))
boxplot(data$human_footprint~phylum_name,data=data, horizontal=TRUE,
        las=2, ylab="",xlab="Human footprint index",col=pal8,frame=F,cex.lab=1.2,cex.axis=1.2,cex=0.8, reverse=T)
dev.off()

###Boxplot of human footprint by COSEWIC assessments for Canada dataset####
bold_can <- filter(data, country == "Canada" & assessment != "") %>% droplevels()
bold_can$assessment<-factor(bold_can$assessment,levels=c("Extinct","Extirpated","Endangered","Threatened","Special Concern","Not At Risk","Data Deficient","Not Available"))

png("../04_figures/20220704_boxplot_footprint_canada.png", width = 4,height = 7, units = "in", res=300)
par(mar=c(4,8.4,0.5,0.32))
boxplot(canada_threat~assessment,data=bold_can, horizontal=TRUE,
        las=2, ylab="",xlab="Canada threat index",col=rev(pal8),frame=F,cex.lab=.2,cex.axis=1.2,cex=0.8, reverse=T)
dev.off()


###### ANOVAS for variations of footprint by different factors ############
###Global data, footprint vs phylum
fit_global<-aov(human_footprint~phylum_name,data=data)
summary(fit_global)

TukeyHSD(fit)
BOLD_global_HSD<-glht(fit_global, linfct=mcp(phylum_name="Tukey"))
cld(BOLD_global_HSD)
bold_g_let<-c("ab", "a", "d","a", "abc", "b","cd")
boxplot(human_footprint~phylum_name, col="#FDBF6F", boxwex = 0.7, cex.axis=.9, ylab="Human footprint index", data=data, ylim=c(0,55), xlab="Taxon")
text(c(1:7), c(rep(54,7)), labels = bold_g_let, cex = 1.2)


###Canadian data, footprint vs COSEWIC
fit_can<-aov(human_footprint~assessment,data=bold_can)
summary(fit_can)

TukeyHSD(fit_can)
BOLD_can_HSD<-glht(fit_can, linfct=mcp(assessment="Tukey"))
cld(BOLD_can_HSD)
bold_can_let<-c("abcd", "abcd", "cd","d", "a", "b","b","c")
boxplot(human_footprint~assessment, col="#FDBF6F", boxwex = 0.7, cex.axis=.9, ylab="Canada's Human Footprint", data=bold_can, ylim=c(0,55), xlab="Taxon")
text(c(1:8), c(rep(54,7)), labels = bold_can_let, cex = 1.2)





####### Variation of human footprint index by taxon, both ANOVAs and Figures

fit_global<-aov(human_footprint~taxon,data=data)
summary(fit_global)
TukeyHSD(fit)
BOLD_global_HSD<-glht(fit_global, linfct=mcp(taxon="Tukey"))
cld(BOLD_global_HSD)
bold_g_let<-c("ab", "a", "d","a", "abc", "b","cd")



fit_can<-aov(human_footprint~taxon,data=bold_can)
summary(fit_can)
TukeyHSD(fit_can)
BOLD_can_HSD<-glht(fit_can, linfct=mcp(taxon="Tukey"))
cld(BOLD_can_HSD)
bold_can_let<-c("c", "a", "d","b", "abcd", "c","abc")


png("../04_figures/20220704_boxplot_footprint_taxon_all.png", width = 8,height = 7, units = "in", res=300)
par(mfrow=c(2,1),mar=c(5,5,0.5,0.32))
boxplot(human_footprint~taxon, col="#FDBF6F", boxwex = 0.7, cex.axis=.9, ylab="Global Human Footprint", data=data, ylim=c(0,55), xlab="Taxon")
text(c(1:7), c(rep(54,7)), labels = bold_g_let, cex = 1.2)
boxplot(human_footprint~taxon, col="#FDBF6F", boxwex = 0.7, cex.axis=.9, ylab="Canada's Human Footprint", data=bold_can, ylim=c(0,55), xlab="Taxon")
text(c(1:7), c(rep(54,7)), labels = bold_can_let, cex = 1.2)
dev.off()


