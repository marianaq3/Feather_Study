# **Feather morphology of Colombian birds**

## **Description of the data and file structure**

We examined 2504 feathers from 1166 Neotropical birds of 444 species, comprising 43 out of the approximately 67 families that can be found in Colombia. Feathers were sampled from birds captured in mist nets at Parque Nacional Natural Tatamá, Parque Nacional Natural Farallones, and Remedios, Antioquia or from museum specimens housed at the Zoological Collections of Universidad Icesi (ICZA), the Museo de Historia Natural de la Universidad de los Andes (ANDES-O), and the bird collection of the Instituto de Investigación de Recursos Biológicos Alexander von Humboldt (IAvH-A).These locations provided a wide variety of precipitation and temperature conditions across elevations ranging from 600 to 2500 masl.

We sampled breast, back and head feathers. Breast feathers were collected from the upper breast region near the pectoralis and supracoracoideus muscles. Back and head feathers were collected from the approximate middle section for each of these body parts. All feathers were flattened and photographed using a Canon T6 camera with a 100 mm Canon macro lens and a Canon macro ring lite MR-14EX (ring flash). 

From each contour feather we measured in millimeters 1) total feather length; 2) length of the plumulaceous downy region; 3) and length of pennaceous distal section, and in squared millimeters, we measured 4) area of the downy region; 5) area of the pennaceous region and 6) total feather area. Measurements were obtained from tip to calamus using the software ImageJ version 2.14.0/1.54f generating straightened images of contour feathers. Using these raw measurements, we calculated for each contour feather 1) the proportion of the downy region, quantified as the length of the downy region in millimeters divided by the total feather length; 2) the proportion of the pennaceous region obtained as the length of the pennaceous region divided by the total feather length; 3) the proportion of the area of the downy region of the feather, quantified as the area of the downy region divided by the total feather area, and 4) the proportion of the area of the pennaceous region of the feather, quantified as the area of the pennaceous region divided by the total feather area. Also, for body mass, we used Ocampo et al. (2021).  

Precipitation and temperature were downloaded from WorldClim version 2.1. Temperature data (bioclimatic variable 12) and precipitation data (bioclimatic variable 1) were downloaded for all locations from the years 1970 to 2000 at a spatial resolution of 2.5 minutes. Mean temperature and precipitation values were calculated using WorldClim.

Data were imported and matched to the coordinates of each feather’s collection location using R packages sp (Bivand et al., 2013; Pebesma & Bivand, 2005), raster (Hijmans, 2020)  and rgdal (Bivand et al., 2021). 

### **Files and variables**

#### **File:** Feather_Models.R

**Description:** R script with all of the codes used to run the phylogenetic least squares models and the graphs developed for the manuscript. 

#### **File:** Tree_2025.tre

**Description:** Phylogenetic tree with the species in the study. This tree was used for the phylogenetic least squares models and to test phylogenetic  signal.

#### **File:** Body_parts_dataset.csv

**Description:** feather measurements for body part specific phylogenetic least squares models. It contains measurements from breast, back and head feathers. 

**Variables**

In_Tree: Species name in phylogenetic tree

Mass: Species average body weight in grams 

FL_Breast: Total breast feather length in millimeters

DL_Breast: Length of downy section of breast feather in millimeters

PDL_Breast: proportion of downy section length of breast feather

PPL_Breast: proportion of pennaceous section length of breast feather

AD_Breast: area of downy section of breast feather in square millimeters

AP_Breast: area of pennaceous section of breast feather in square millimeters

AT_Breast: total area of breast feather in square millimeters

PAD_Breast: proportion of downy section area of breast feather

PAP_Breast: proportion of pennaceous section area of breast feather

FL_Back: Total back feather length in millimeters

DL_Back: Length of downy section of back feather in millimeters

PDL_Back: proportion of downy section length of back feather

PPL_Back: proportion of pennaceous section length of back feather

AD_Back: area of downy section of back feather in square millimeters

AP_Back: area of pennaceous section of back feather in square millimeters

AT_Back: total area of back feather in square millimeters

PAD_Back: proportion of downy section area of back feather

PAP_Back: proportion of pennaceous section area of back feather

FL_Head: Total head feather length in millimeters

DL_Head: Length of downy section of head feather in millimeters

PDL_Head: proportion of downy section length of head feather

PPL_Head: proportion of pennaceous section length of head feather

AD_Head: area of downy section of head feather in square millimeters

AP_Head: area of pennaceous section of head feather in square millimeters

AT_Head: total area of head feather in square millimeters

PAD_Head: proportion of downy section area of head feather

PAP_Head: proportion of pennaceous section area of head feather

Temperature: temperature in degrees celsius

Rainfall: millimeters of rain

#### **File:** General_model_dataset.csv

**Description:** feather measurements for general body feather phylogenetic least squares model. 

**Variables**

In_Tree: Species name in phylogenetic tree

Mass: Species average body weight in grams

FL: Total feather length in millimeters

DL: Length of downy section in millimeters

PL: Length of pennaceous section in millimeters

PDL: Proportion of downy section length

PPL: Proportion of pennaceous section length

AD: Area of downy section in square millimeters

AP: Area of pennaceous section in square millimeters

AT: Total feather area in square millimeters

PAD: Proportion of downy section area

PAP: Proportion of pennaceous section area

Temperature: temperature in degrees celsius

Rainfall: millimeters of rain

## Code/software

R and RStudio are needed to view and run the R script used for the data analysis and plot generation. All packages and workflows are described on the script. 
