## Packages needed

library(lme4)
library(MuMIn)
library(nlme)
library(phytools)
library(ape)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(effects)
library(gridExtra)
library(sf)
library(mapview)
library(webshot2)
library(grid)
library(geiger)
library(devtools)
library(patchwork)
library(RColorBrewer)
library(raster)
library(gstat)
library(sp)
library(rnaturalearth)
library(ggtree)
library(rsvg)
library(cowplot)


## Data upload and type transformation
setwd("~/Desktop/")
Feathers_All <- read.csv("General_Model.csv", 
                         header = T)
str(Feathers_All)

transform_to_numeric <- function(df, colnames_to_convert) {
  df[colnames_to_convert] <- lapply(df[colnames_to_convert], as.numeric)
  return(df)
}

Feathers_All <- transform_to_numeric(Feathers_All,c("Rainfall","AP",
                                                    "AT","PAD",
                                                    "PAP"))
str(Feathers_All)

sub.Feathers_All	<- Feathers_All[complete.cases(Feathers_All),]

## Scale data so that all numbers are in the same order of magnitude
Feathers_All_sc <- data.frame(Order=Feathers_All$Order,
                              Family=Feathers_All$Family,
                              Genus=Feathers_All$Genus,
                              Species=Feathers_All$Current_Identification,
                              In_Tree=Feathers_All$In_Tree,
                              Mass=scale(Feathers_All$Mass),
                              FL=scale(Feathers_All$FL),
                              DL=scale(Feathers_All$DL),
                              PL=scale(Feathers_All$PL),
                              PDL=scale(Feathers_All$PDL),
                              PPL=scale(Feathers_All$PPL), 
                              AD=scale(Feathers_All$AD),
                              AP=scale(Feathers_All$AP),
                              AT=scale(Feathers_All$AT),
                              PAD=scale(Feathers_All$PAD),
                              PAP=scale(Feathers_All$PAP),
                              Temperature=scale(Feathers_All$Temperature),
                              Rainfall=scale(Feathers_All$Rainfall))
str(Feathers_All_sc)

sub.Feathers_All_sc	<-	Feathers_All_sc[complete.cases(Feathers_All_sc),]

## Phylogenetic test
# Number of significant digits
significant_digits = 3

# Calculate mean by group for specified columns
mean_by_sp <- sub.Feathers_All_sc %>%
  group_by(In_Tree) %>%
  summarise(across(c(Mass,FL,DL,PL,PDL,PPL,
                     AD,AP,AT,PAD,PAP,Temperature,
                     Rainfall), ~ signif(mean(.), significant_digits)))

str(mean_by_sp)


write.table(mean_by_sp,"~/Desktop/mean_by_sp.csv",sep=",",row.names=FALSE)

#Import tree
tree1<- read.nexus(paste(getwd(),"Tree_2025.tre",sep="/"))
Ntip(tree1)
spp.vec	<-	levels(as.factor(mean_by_sp$In_Tree))
spp.vec	<-	gsub(" ","_",mean_by_sp$In_Tree)

# Check that the tips of th tree and names in data set are the same
sum(!spp.vec%in%tree1$tip.label)
spp.vec[!spp.vec%in%tree1$tip.label]
spp.vec[!tree1$tip.label%in%spp.vec]
match(tree1$tip.label,spp.vec)
match(spp.vec,tree1$tip.label)
length(spp.vec)

# Organize data in phylogenetic order
feathers.phy	<-	mean_by_sp[match(tree1$tip.label,spp.vec),]
rownames(feathers.phy)	<-	gsub(x = feathers.phy$In_Tree, pattern = " ", 
                               replacement = "_")

is.binary.tree(tree1)
is.ultrametric(tree1)

write.table(feathers.phy,"~/Desktop/feathers.phy.csv",sep=",",
            row.names=FALSE)

str(feathers.phy)

#Phylogenetic signal evaluation
phylosig(tree1, feathers.phy$PAD, method = "lambda", test = TRUE)


## PGLS
## General model
phy_mod_PAD	<-	gls(PAD~Mass+Temperature+Rainfall+AT+
                     Temperature*Rainfall,
                   data=feathers.phy,
                   correlation=corBrownian(value=1,phy=tree1),
                   method="ML")

phy_mod_PAD

# Model selection
mod_sel_PAD	<-	dredge(phy_mod_PAD,rank="BIC")
mod_sel_PAD

## Models for each body part
# Import data
Feathers_BHB <- read.csv("Feathers_BHB.csv", 
                         header = T)
str(Feathers_BHB)

transform_to_numeric <- function(df, colnames_to_convert) {
  df[colnames_to_convert] <- lapply(df[colnames_to_convert], as.numeric)
  return(df)
}

Feathers_BHB <- transform_to_numeric(Feathers_BHB,c("Rainfall","AP_Head",
                                                    "AT_Head","PAD_Head",
                                                    "PAP_Head"))
str(Feathers_BHB)

sub.Feathers_BHB	<- Feathers_BHB[complete.cases(Feathers_BHB),]

# Scale data frame
Feathers_BHB_sc <- data.frame(Order=Feathers_BHB$Order,
                              Family=Feathers_BHB$Family,
                              Genus=Feathers_BHB$Genus,
                              Species=Feathers_BHB$Current_Identification,
                              In_Tree=Feathers_BHB$In_Tree,
                              Mass=scale(Feathers_BHB$Mass),
                              FL_Breast=scale(Feathers_BHB$FL_Breast),
                              DL_Breast=scale(Feathers_BHB$DL_Breast),
                              PDL_Breast=scale(Feathers_BHB$PDL_Breast),
                              PPL_Breast=scale(Feathers_BHB$PPL_Breast), 
                              AD_Breast=scale(Feathers_BHB$AD_Breast),
                              AP_Breast=scale(Feathers_BHB$AP_Breast),
                              AT_Breast=scale(Feathers_BHB$AT_Breast),
                              PAD_Breast=scale(Feathers_BHB$PAD_Breast),
                              PAP_Breast =scale(Feathers_BHB$PAP_Breast),
                              FL_Back=scale(Feathers_BHB$FL_Back),
                              DL_Back=scale(Feathers_BHB$DL_Back),
                              PDL_Back=scale(Feathers_BHB$PDL_Back),
                              PPL_Back=scale(Feathers_BHB$PPL_Back), 
                              AD_Back=scale(Feathers_BHB$AD_Back),
                              AP_Back=scale(Feathers_BHB$AP_Back),
                              AT_Back=scale(Feathers_BHB$AT_Back),
                              PAD_Back=scale(Feathers_BHB$PAD_Back),
                              PAP_Back=scale(Feathers_BHB$PAP_Back),
                              FL_Head=scale(Feathers_BHB$FL_Head),
                              DL_Head=scale(Feathers_BHB$DL_Head),
                              PDL_Head=scale(Feathers_BHB$PDL_Head),
                              PPL_Head=scale(Feathers_BHB$PPL_Head), 
                              AD_Head=scale(Feathers_BHB$AD_Head),
                              AP_Head=scale(Feathers_BHB$AP_Head),
                              AT_Head=scale(Feathers_BHB$AT_Head),
                              PAD_Head=scale(Feathers_BHB$PAD_Head),
                              PAP_Head=scale(Feathers_BHB$PAP_Head),
                              Temperature=scale(Feathers_BHB$Temperature),
                              Rainfall=scale(Feathers_BHB$Rainfall))
str(Feathers_BHB_sc)

sub.Feathers_BHB_sc	<-	Feathers_BHB_sc[complete.cases(Feathers_BHB_sc),]

# Number of significant digits
significant_digits = 3

# Calculate mean by group for specified columns
mean_by_sp_body <- sub.Feathers_BHB_sc %>%
  group_by(In_Tree) %>%
  summarise(across(c(Mass,FL_Breast,DL_Breast,PDL_Breast,PPL_Breast,
                     AD_Breast,AP_Breast,AT_Breast,PAD_Breast,PAP_Breast,FL_Back,
                     DL_Back,PDL_Back,PPL_Back,AD_Back,AP_Back,AT_Back,
                     PAD_Back,PAP_Back,FL_Head,DL_Head,PDL_Head,PPL_Head,
                     AD_Head,AP_Head,AT_Head,PAD_Head,PAP_Head,Temperature,
                     Rainfall), ~ signif(mean(.), significant_digits)))

str(mean_by_sp_body)


write.table(mean_by_sp_body,"~/Desktop/mean_by_sp_body",sep=",",row.names=FALSE)

# Match data with tree tips
Ntip(tree1)
spp.vec.body	<-	levels(as.factor(mean_by_sp_body$In_Tree))
spp.vec.body	<-	gsub(" ","_",mean_by_sp_body$In_Tree)

# Check that the tips of th tree and names in data set are the same
sum(!spp.vec.body%in%tree1$tip.label)
spp.vec.body[!spp.vec.body%in%tree1$tip.label]
spp.vec.body[!tree1$tip.label%in%spp.vec.body]
match(tree1$tip.label,spp.vec.body)
match(spp.vec.body,tree1$tip.label)
length(spp.vec.body)

# Organize data in phylogenetic order
feathers.phy_body	<-	mean_by_sp_body[match(tree1$tip.label,spp.vec.body),]

sub.feathers.phy.body	<-	feathers.phy_body[complete.cases(feathers.phy_body),]

rownames(sub.feathers.phy.body)	<-	gsub(x = sub.feathers.phy.body$In_Tree, pattern = " ", 
                                        replacement = "_")

write.table(sub.feathers.phy.body,"~/Desktop/sub.feathers.phy.body.csv",sep=",",
            row.names=FALSE)

str(sub.feathers.phy.body)

# Phylogenetic signal for breast
phylosig(tree1, feathers.phy_body$PAD_Breast, method = "lambda", test = TRUE)

# Phylogenetic signal for back
phylosig(tree1, feathers.phy_body$PAD_Back, method = "lambda", test = TRUE)

# Phylogenetic signal for head
phylosig(tree1, feathers.phy_body$PAD_Head, method = "lambda", test = TRUE)

# PGLS Breast
phy_mod_PAD_Br	<-	gls(PAD_Breast~Mass+Temperature+Rainfall+AT_Breast+
                        Temperature*Rainfall,
                      data=sub.feathers.phy.body,
                      correlation=corBrownian(value=1,phy=tree1),
                      method="ML")

phy_mod_PAD_Br

# Model selection Breast
mod_sel_PAD_Br	<-	dredge(phy_mod_PAD_Br,rank="BIC")
mod_sel_PAD_Br

# PGLS Back
phy_mod_PAD_Ba	<-	gls(PAD_Back~Mass+Temperature+Rainfall+AT_Back+
                        Temperature*Rainfall,
                      data=sub.feathers.phy.body,
                      correlation=corBrownian(value=1,phy=tree1),
                      method="ML")

phy_mod_PAD_Ba

# Model selection Back
mod_sel_PAD_Ba	<-	dredge(phy_mod_PAD_Ba,rank="BIC")
mod_sel_PAD_Ba

# PGLS Head
phy_mod_PAD_He	<-	gls(PAD_Head~Mass+Temperature+Rainfall+AT_Head+
                        Temperature*Rainfall,
                      data=sub.feathers.phy.body,
                      correlation=corBrownian(value=1,phy=tree1),
                      method="ML")

phy_mod_PAD_He

# Model selection Head
mod_sel_PAD_He	<-	dredge(phy_mod_PAD_He,rank="BIC")
mod_sel_PAD_He

## Model graphs
## BREAST
# Effect of Mass
mass_effect_PGLS <- effect("Mass", phy_mod_PAD_Br)
mass_br_PGLS <- plot(mass_effect_PGLS, main="",
                     xlab = "Mass (g)", ylab = "Proportion of downy area")
mass_br_PGLS

# Effect of Temperature
temp_effect_PGLS <- effect("Temperature", phy_mod_PAD_Br)
temp_br_PGLS <- plot(temp_effect_PGLS, main="",
                     xlab = "Temperature (ºC)", ylab = "Proportion of downy area")
temp_br_PGLS

# Effect of Rainfall
rain_effect_PGLS <- effect("Rainfall", phy_mod_PAD_Br)
rain_br_PGLS <- plot(rain_effect_PGLS, main="",
                     xlab = "Rainfall (mm)", ylab = "Proportion of downy area")
rain_br_PGLS

# Effect of AT_Breast
atbreast_effect_PGLS <- effect("AT_Breast", phy_mod_PAD_Br)
at_br_PGLS <- plot(atbreast_effect_PGLS, main="", 
                   xlab = "Total area of breast feather (mm2)",
                   ylab = "Proportion of downy area")
at_br_PGLS

# Interaction Effect of Rainfall and Temperature
interaction_effect_PGLS <- effect("Temperature:Rainfall", phy_mod_PAD_Br,
                                  xlab = "Temperature:Rainfall", ylab = "Proportion of downy area")
temp_rain_br_PGLS <- plot(interaction_effect_PGLS, main="", ylab = "Proportion of downy area")
temp_rain_br_PGLS

par(cex.axis=20, cex.lab=20, cex.main=20)
breast_PGLS <- grid.arrange(mass_br_PGLS, temp_br_PGLS,at_br_PGLS, ncol=3)

# BACK
# Effect of Mass
mass_effect_PGLS <- effect("Mass", phy_mod_PAD_Ba)
mass_ba_PGLS <- plot(mass_effect_PGLS, main="",
                     xlab = "Mass (g)", ylab = "Proportion of downy area")
mass_ba_PGLS

# Effect of Temperature
temp_effect_PGLS <- effect("Temperature", phy_mod_PAD_Ba)
temp_ba_PGLS <- plot(temp_effect_PGLS, main="",
                     xlab = "Temperature (ºC)", ylab = "Proportion of downy area")
temp_ba_PGLS

# Effect of Rainfall
rain_effect_PGLS <- effect("Rainfall", phy_mod_PAD_Ba)
rain_ba_PGLS <- plot(rain_effect_PGLS, main="",
                     xlab = "Rainfall (mm)", ylab = "Proportion of downy area")
rain_ba_PGLS

# Effect of AT_Back
atback_effect_PGLS <- effect("AT_Back", phy_mod_PAD_Ba)
at_ba_PGLS <- plot(atback_effect_PGLS, main="", 
                   xlab = "Total area of back feather (mm2)",
                   ylab = "Proportion of downy area")
at_ba_PGLS

# Interaction Effect of Rainfall and Temperature
interaction_effect_PGLS <- effect("Temperature:Rainfall", phy_mod_PAD_Ba,
                                  xlab = "Temperature:Rainfall", 
                                  ylab = "Proportion of downy area")
temp_rain_ba_PGLS <- plot(interaction_effect_PGLS, 
                          main="",ylab = "Proportion of downy area")
temp_rain_ba_PGLS

par(cex.axis=20, cex.lab=20, cex.main=20)
back_PGLS <- grid.arrange(temp_ba_PGLS,rain_ba_PGLS,temp_rain_ba_PGLS, ncol=2)

layout_matrix <- rbind(c(1, 2),
                       c(3, 3))

back_PGLS2 <- grid.arrange(
  temp_ba_PGLS, 
  rain_ba_PGLS, 
  temp_rain_ba_PGLS, 
  layout_matrix = layout_matrix,
  heights = c(1, 1.4)  # give more vertical room to bottom row
)

back_PGLS2

## HEAD
# Effect of Mass
mass_effect_PGLS <- effect("Mass", phy_mod_PAD_He)
mass_he_PGLS <- plot(mass_effect_PGLS, main="",
                     xlab = "Mass (g)", ylab = "Proportion of downy area")
mass_he_PGLS

# Effect of Temperature
temp_effect_PGLS <- effect("Temperature", phy_mod_PAD_He)
temp_he_PGLS <- plot(temp_effect_PGLS, main="",
                     xlab = "Temperature (ºC)", ylab = "Proportion of downy area")
temp_he_PGLS

# Effect of Rainfall
rain_effect_PGLS <- effect("Rainfall", phy_mod_PAD_He)
rain_he_PGLS <- plot(rain_effect_PGLS, main="",
                     xlab = "Rainfall (mm)", ylab = "Proportion of downy area")
rain_he_PGLS

# Effect of AT_Head
athead_effect_PGLS <- effect("AT_Head", phy_mod_PAD_He)
at_he_PGLS <- plot(athead_effect_PGLS, main="", 
                   xlab = "Total area of head feather (mm2)",
                   ylab = "Proportion of downy area")
at_he_PGLS

# Interaction Effect of Rainfall and Temperature
interaction_effect_PGLS <- effect("Temperature:Rainfall", phy_mod_PAD_He,
                                  xlab = "Temperature:Rainfall", 
                                  ylab = "Proportion of downy area")
temp_rain_he_PGLS <- plot(interaction_effect_PGLS, 
                          main="",ylab = "Proportion of downy area")
temp_rain_he_PGLS

par(cex.axis=20, cex.lab=20, cex.main=20)
head_PGLS <- grid.arrange(mass_he_PGLS, temp_he_PGLS,at_he_PGLS, ncol=3)

## ALLOMETRIC AND OTHER PLOTS
log.rawmass <- log(sub.Feathers_All$Mass)
log.rawlength <- log(sub.Feathers_All$FL)
sub.Feathers_All$Log_Length <- log.rawlength
sub.Feathers_All$Log_Mass <- log.rawmass
sub.Feathers_All$rel_mass <- sub.Feathers_All$FL/sub.Feathers_All$Mass

# Feather length - general
FL_Mass <- ggscatter(sub.Feathers_All, x= "Log_Mass", y= "Log_Length",
                     add = "reg.line", conf.int = T,
                     cor.coef = T, cor.method = "pearson", cor.coef.size = 5, 
                     cor.coef.coord = c(1,5),
                     xlab = "Avian weight (g)", 
                     ylab = "Feather length (mm)", 
                     color = "Temperature") + 
  theme(text = element_text(size=16), legend.position = "top")+  
  scale_color_gradient(low = "blue", high = "red")+
  labs(colour="Temperature ºC") +
  annotate(family="serif",geom = "text", x = 1, y = 1.55, 
           label = "A", size = 8, vjust = 1.55,col="black",fontface="bold")
FL_Mass

# Breast
log.rawmass <- log(Feathers_BHB$Mass)
log.rawlength <- log(Feathers_BHB$FL_Breast)
Feathers_BHB$Log_Length_BR <- log.rawlength
Feathers_BHB$Log_Mass <- log.rawmass
Feathers_BHB$rel_mass <- Feathers_BHB$FL_Breast/Feathers_BHB$Mass

Br_FL_Mass <- ggscatter(Feathers_BHB, x= "Log_Mass", y= "Log_Length_BR",
                        add = "reg.line", conf.int = T,
                        cor.coef = T, cor.method = "pearson", cor.coef.size = 5, 
                        cor.coef.coord = c(1,5),
                        xlab = "", 
                        ylab = "Feather length (mm): Breast", 
                        color = "Temperature") + 
  theme(text = element_text(size=16), legend.position = "top")+  
  scale_color_gradient(low = "blue", high = "red")+
  labs(colour="Temperature ºC") +
  annotate(family="serif",geom = "text", x = 1, y = 1.55, 
           label = "A", size = 8, vjust = 1.55,col="black",fontface="bold")
Br_FL_Mass

# Back
log.rawmass1 <- log(sub.Feathers_BHB$Mass)
log.rawlength1 <- log(sub.Feathers_BHB$FL_Back)
sub.Feathers_BHB$Log_Length_Ba <- log.rawlength1
sub.Feathers_BHB$Log_Mass1 <- log.rawmass1
sub.Feathers_BHB$rel_mass1 <- sub.Feathers_BHB$FL_Back/sub.Feathers_BHB$Mass

Ba_FL_Mass <- ggscatter(sub.Feathers_BHB, x= "Log_Mass1", y= "Log_Length_Ba",
                        add = "reg.line", conf.int = T,
                        cor.coef = T, cor.method = "pearson", cor.coef.size = 5, 
                        cor.coef.coord = c(1,5),
                        xlab = "", 
                        ylab = "Feather length (mm): Back", 
                        color = "Temperature") + 
  theme(text = element_text(size=16), legend.position = "top")+  
  scale_color_gradient(low = "blue", high = "red")+
  labs(colour="Temperature ºC") +
  annotate(family="serif",geom = "text", x = 1, y = 1.55, 
           label = "B", size = 8, vjust = 1.55,col="black",fontface="bold")
Ba_FL_Mass

# Head
log.rawmass2<- log(sub.Feathers_BHB$Mass)
log.rawlength2 <- log(sub.Feathers_BHB$FL_Head)
sub.Feathers_BHB$Log_Length_He <- log.rawlength2
sub.Feathers_BHB$Log_Mass2 <- log.rawmass2
sub.Feathers_BHB$rel_mass2 <- sub.Feathers_BHB$FL_Head/sub.Feathers_BHB$Mass

He_FL_Mass <- ggscatter(sub.Feathers_BHB, x= "Log_Mass2", y= "Log_Length_He",
                        add = "reg.line", conf.int = T,
                        cor.coef = T, cor.method = "pearson", cor.coef.size = 5, 
                        cor.coef.coord = c(1,5),
                        xlab = "", 
                        ylab = "Feather length (mm): Head", 
                        color = "Temperature") + 
  theme(text = element_text(size=16), legend.position = "top")+  
  scale_color_gradient(low = "blue", high = "red")+
  labs(colour="Temperature ºC") +
  annotate(family="serif",geom = "text", x = 1, y = 1.55, 
           label = "C", size = 8, vjust = 1.55,col="black",fontface="bold")
He_FL_Mass

# Merge allometric plots
FL_Mass_All <- ggarrange(Br_FL_Mass,Ba_FL_Mass,
                         He_FL_Mass,
                         ncol = 1, nrow=3, font.label=list(color="black",size=15,
                                                           face="plain"),
                         common.legend = TRUE, hjust = -3.5, vjust = 1)

FL_Mass_All
FL_Mass_All1 <- annotate_figure(FL_Mass_All, bottom = textGrob("Logarithm of the average body weight (g)", 
                                                               just = "centre"))
FL_Mass_All1

## Rainfall and Temperature plots
# Rainfall
sub.Feathers_All	<- Feathers_All[complete.cases(Feathers_All),]

# Define rainfall groups
sub.Feathers_All$Rainfall_Group <- cut(sub.Feathers_All$Rainfall,
                                       breaks = c(700, 1500, 2500, 3500, 3999, Inf),
                                       labels = c("700-1500 mm", "1501-2500 mm", "2501-3500 mm", "3501-3999 mm", "4000+ mm"),
                                       include.lowest = TRUE)

rainfall.plot2 <- ggboxplot(sub.Feathers_All, x = "Rainfall_Group", y = "PDL",
                            xlab = "Rainfall Group (mm)", 
                            ylab = "Proportion of Downy", 
                            fill = "Rainfall_Group") + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +  # Add jittered data points
  scale_fill_manual(values = c("lightskyblue1", "skyblue3", "royalblue3", "mediumblue", "darkblue")) +
  ylim(0.2, 0.75) +
  theme_minimal() +
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

rainfall.plot2

# Temperature
# Define temperature groups
sub.Feathers_All$Temperature_Group <- cut(sub.Feathers_All$Temperature,
                                          breaks = c(9, 15, 22, 28),
                                          labels = c("9-15°C", "16-22°C", "23-28°C"),
                                          include.lowest = TRUE)

sub.Feathers_All <- sub.Feathers_All[!is.na(sub.Feathers_All$Temperature_Group), ]

Temp.plot <- ggboxplot(sub.Feathers_All, x = "Temperature_Group", y = "PDL",
                       xlab = "Temperature Group (ºC)", 
                       ylab = "Proportion of Downy", 
                       fill = "Temperature_Group") + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +  # Add jittered data points
  scale_fill_manual(values =  c("lightblue", "orange", "red")) +
  ylim(0.2, 0.75) +
  theme_minimal() +
  theme(text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

Temp.plot

# One plot for both environmental conditions
Temp.plot <- Temp.plot + theme(legend.position = "none")
rainfall.plot2 <- rainfall.plot2 + theme(legend.position = "none")

Temp_Rain <- ggarrange(Temp.plot,rainfall.plot2,
                       ncol = 2, nrow=1, font.label=list(color="black",size=15,
                                                         face="plain"),
                       hjust = -3.5, vjust = 1)

Temp_Rain

# Temperature and Rainfall maps
# Load Colombia boundary and project
colombia <- ne_countries(scale = "medium", country = "Colombia", returnclass = "sf")
colombia_proj <- st_transform(colombia, crs = 3116)
colombia_sp <- as(colombia_proj, "Spatial")

# Create a raster grid from Colombia's extent
r_template <- raster(extent(colombia_sp), res = 10000, crs = proj4string(colombia_sp))
r_masked <- mask(r_template, colombia_sp)  # Raster masked to Colombia shape

# Rasterize Colombia shape
r_colombia <- rasterize(colombia_sp, r_template, field = 1)

# make sure it's not all NA
if (all(is.na(values(r_colombia)))) {
  stop("Colombia rasterization failed — check CRS alignment.")
}

# Use SpatialPixels from the valid raster cells
grid_coords <- as(r_colombia, "SpatialPixelsDataFrame")
idw_temp <- idw(Temperature ~ 1, data_sp, newdata = grid_coords)
temp_raster <- raster(idw_temp)

temp_clipped <- mask(temp_raster, colombia_sp)

temp_df <- as.data.frame(rasterToPoints(temp_clipped))
names(temp_df) <- c("x", "y", "Temperature")

# IDW interpolation for Rainfall using same grid
idw_rain <- idw(Rainfall ~ 1, data_sp, newdata = grid_coords)
rain_raster <- raster(idw_rain)

# Mask to Colombia boundary
rain_clipped <- mask(rain_raster, colombia_sp)

# Convert to data frame for plotting
rain_df <- as.data.frame(rasterToPoints(rain_clipped))
names(rain_df) <- c("x", "y", "Rainfall")

# Temperature Map
temp_map <- ggplot() +
  geom_raster(data = temp_df, aes(x = x, y = y, fill = Temperature)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), name = "Temperature (°C)") +
  geom_sf(data = data_proj, color = "black", size = 0.8) +
  geom_sf(data = colombia_proj, fill = NA, color = "black", size = 0.6) +
  coord_sf(crs = st_crs(3116)) +
  theme_minimal() + labs(x = "Longitude (°W)", y = "Latitude (°N)")


# Rainfall Map
rain_map <- ggplot() +
  geom_raster(data = rain_df, aes(x = x, y = y, fill = Rainfall)) +
  scale_fill_gradientn(colors = c("#d0d1e6", "#74a9cf", "#0570b0", "#023858"), name = "Rainfall (mm)")+
  geom_sf(data = data_proj, color = "black", size = 0.8) +
  geom_sf(data = colombia_proj, fill = NA, color = "black", size = 0.6) +
  coord_sf(crs = st_crs(3116)) +
  theme_minimal() + labs(x = "Longitude (°W)", y = "Latitude (°N)")

# Combined map
combined_map <- temp_map + rain_map +
  plot_layout(ncol = 2) 

combined_map