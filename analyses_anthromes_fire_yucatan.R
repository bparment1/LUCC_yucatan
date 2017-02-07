#################################    FIRE YUCATAN PAPER  #######################################
#############################         LAND COVER CHANGE          #######################################
#This script reads data extracted from MODIS FIRE time series                             
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial multinomial logit model with explantory variables.                                                     #
#The link between fire as a tool for anthropogenic changes is also explore with the anthrome dataset.
#We used a non negative binomial model with Tukey to examine differences between anthromes classes.
#
#AUTHOR: Benoit Parmentier                                                                #
#DATE CREATED: 02/06/2016 
#DATE MODIFIED: 02/07/2017
#Version: 1
#PROJECT: Land cover Change Yucatan with Marco Millones 
#   
#TODO:
#COMMENTS:
#
#COMMIT: clean up of code and check of results for negative binomial model

#################################################################################################

#################################################################################################

###Loading R library and packages                                                      
library(gtools)                          # loading some useful tools 
library(mgcv)                            # GAM package by Simon Wood
library(sp)                              # Spatial pacakge with class definition by Bivand et al.
library(spdep)                           # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                           # GDAL wrapper for R, spatial utilities
library(gstat)                           # Kriging and co-kriging by Pebesma et al.
library(fields)                          # NCAR Spatial Interpolation methods such as kriging, splines
library(raster)                          # Hijmans et al. package for raster processing
library(foreign)                         # Library for format exchange (e.g. dbf,spss,sas etc.)
library(gdata)                           # various tools with xls reading
library(xts)                             # basic package for time series analysis
library(zoo)                             # basic package for time series analysis
#library(forecast)                       # package containing ARIMA procedures
library(rasterVis)                       # plotting raster
library(nnet)                            # Contains multinom and neural net functions
library(ggplot2)                         # plotting package
library(reshape2)                        # data wrangling
library(mlogit)                          # maximum liklihood estimation and multinomial model
library(parallel)                        # parralel programming and multi cores
library(plyr)
library(rgeos)                           # topology and vector spatial queries and operations
library(afex)                            # functions related to ANOVA
library(car)
library(MASS)                            # contains negative binomial model
library(multcomp)                        # contains Tukey comparison

###### Functions used in this script

functions_analyses_script <- "analyses_fire_yucatan_functions_05282016.R" #PARAM 1
functions_analyses_anthromes_script <- "analyses_anthromes_fire_yucatan_functions_02072017.R" #PARAM 1

script_path <- "/home/bparmentier/Google Drive/FireYuca_2016/R_scripts" #path to script #PARAM 2

source(file.path(script_path,functions_analyses_script)) #source all functions used in this script 1.
source(file.path(script_path,functions_analyses_anthromes_script)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir <- "/home/bparmentier/Google Drive/FireYuca_2016/"
out_dir <- "/home/bparmentier/Google Drive/FireYuca_2016/outputs"

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"yucatan_CI_analyses_anthromes_fire_02012017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
id_name <- "yucatan_fire_pointid" #Column with the reference point id

state_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/IN_QGIS/State_dis_from_muni.shp"
#data_Hansen_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/Hansen_fire.xlsx" #contains the whole dataset
data_Hansen_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/Hansen_fire.csv" #contains the whole dataset

#data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/yucatan_fire.xlsx" #contains the whole dataset
data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/datasets/Firedata_05182016.txt" #contains the whole dataset

#data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/old_data/000_GYR_FIRENDVI_2000-9.txt"
#data_CI_02062016.txt
coord_names <- c("POINT_X","POINT_Y") #PARAM 11
zonal_var_name <- "state" #name of the variable to use to run the model by zone, Yucatan state here
y_var_name <- "fpnfpch"
ref_var_name <- ""
num_cores <- 3
#disable scientific notation
options(scipen = 999)

### Anthropogenic biomes/anthromes inputs:
global_anthromes_filename <- "/home/bparmentier/Google Drive/FireYuca_2016/datasets/gl-anthrome-geotif/gl_anthrome.tif"
#anthromes_legend_filename <- "/home/bparmentier/Google Drive/FireYuca_2016/datasets/gl-anthrome-geotif/anthropogenic_biome_legend.csv"
anthromes_legend_filename <- "/home/bparmentier/Google Drive/FireYuca_2016/datasets/gl-anthrome-geotif/anthromes_legend.csv"
anthromes_yucatan_filename <- "/home/bparmentier/Google Drive/FireYuca_2016/datasets/yucatan_window_anthromes.tif"

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

data_df <-read.table(data_CI_fname,stringsAsFactors=F,header=T,
                     sep=",")
#data_Hansen <-read.table(data_CI_fname,stringsAsFactors=F,header=T,sep=",")

#Create the count variable
fire_modis_col <- c("FIRE_2000","FIRE_2001","FIRE_2002","FIRE_2003","FIRE_2004","FIRE_2005","FIRE_2006","FIRE_2007")
data_df$FIRE_freq <- data_df$FIRE_2000 + data_df$FIRE_2001 + data_df$FIRE_2002 + data_df$FIRE_2003 + data_df$FIRE_2004 + data_df$FIRE_2005 + data_df$FIRE_2006 + data_df$FIRE_2007 
data_df$FIRE_intensity <- data_df$FIRE_freq/8
data_df$FIRE_bool <- data_df$FIRE_freq > 0
data_df$FIRE_bool <- as.numeric(data_df$FIRE_bool)

data_df_spdf <- data_df
coordinates(data_df_spdf) <- coord_names
filename<-sub(extension(basename(state_fname)),"",basename(state_fname))       #Removing path and the extension from file name.
state_outline <- readOGR(dsn=dirname(state_fname), filename)

proj4string(data_df_spdf) <- proj4string(state_outline)
l_poly <- over(data_df_spdf,state_outline) #points in polygon operation
l_poly[[id_name]] <- data_df_spdf[[id_name]]
data_df_spdf<- merge(data_df_spdf,l_poly,by=id_name)
#data_df_spdf <- cbind(data_df_spdf,l_poly) #join back the data
#coordinates(data_df_spdf) <- coord_names
#unique(data_df_spdf$FIRST_NOM_)
data_df_spdf$state <- as.character(data_df_spdf$FIRST_NOM_)
table(data_df_spdf$state)
sum(table(data_df_spdf$state))
nrow(data_df_spdf)

### Quick check of plotting variables
spplot(data_df_spdf,"cattledensity")
spplot(data_df_spdf,"FIRE_pre07")
spplot(data_df_spdf,"dist_roads")
spplot(data_df_spdf,"state")

hist(data_df[[y_var_name]])
table(data_df[[y_var_name]])
barplot(table(data_df[[y_var_name]]))

#### Part II: compare anthromes and fire data ####

#Run for the overall state
out_suffix_s <- paste("overall_",out_suffix,sep="")
df_anthromes_legend <- read.table(anthromes_legend_filename,sep=",",header=T)

#### Crop to the region of interest

#ref_e <- extent(reg_ref_rast) #extract extent from raster object
ref_e <- extent(data_df_spdf) #extract extent from raster object
reg_outline_poly <- as(ref_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
reg_outline_poly <- as(reg_outline_poly, "SpatialPolygonsDataFrame") #promote to spdf
proj4string(reg_outline_poly) <- proj4string(data_df_spdf) #Assign projection to spdf

r_gl_anthromes <- raster(global_anthromes_filename)

reg_outline_poly_WGS84 <- spTransform(reg_outline_poly,CRS(projection(r_gl_anthromes)))
r_test <- crop(r_gl_anthromes,reg_outline_poly_WGS84)
infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="") #name of newly crated shapefile with the extent
writeOGR(reg_outline_poly,dsn= out_dir,layer= sub(".shp","",infile_reg_outline), 
         driver="ESRI Shapefile",overwrite_layer="TRUE")

data_spdf_CRS_WGS84 <- spTransform(data_df_spdf,CRS_WGS84)
#writeOGR(data_spdf_CRS_WGS84,"data_spdf_CRS_WGS84")

outfile<-paste("data_spdf_CRS_WGS84_",out_suffix,sep="")
writeOGR(data_spdf_CRS_WGS84 ,dsn= ".",layer= outfile, driver="ESRI Shapefile",overwrite_layer=TRUE)

#####################
###Anthromes matching the Yucatan peninsula previously cropped in QGIS

r_anthromes_yucatan <- raster(anthromes_yucatan_filename)

#rasterize using sum?
r_fire_freq <- rasterize(data_spdf_CRS_WGS84,y=r_anthromes_yucatan,field="FIRE_freq",fun=sum)

plot(table(as.vector(r_fire_freq)),type="h",
     ylab="Frequency of Fire Freq values in 10km2 pixels",xlab="Fire freq value")

plot(table(as.vector(r_fire_freq)),type="h",
     ylab="Frequency of Fire Freq values in 10km2 pixels",xlab="Fire freq value",
     xlim=c(0,500))
table(as.vector(r_fire_freq))

#df_fire_anthromes <- extract(r_fire_freq,r_anthromes_yucatan)
#tb <- crosstab(r_test,r_anthromes_yucatan)

cat_names <- df_anthromes_legend$LABEL

nb_col <- length(cat_names)
col_anth<-rainbow(nb_col)

plot(r_anthromes_yucatan,col=col_anth,legend=FALSE,axes="FALSE")
legend("topleft",legend=cat_names,title="Anthromes",
       pt.cex=1.1,cex=1.1,fill=col_anth,bty="n")
#scale_position<-c(450000, 600000)
#arrow_position<-c(900000, 600000)
#dev.off()

plot(r_fire_freq,main="Fire frequency count over 2001-2008")

### Reclass: into fewer groups
r_anthromes_yucatan <- raster(r_anthromes_yucatan_filename)
freq(r_anthromes_yucatan)
df_freq <- as.data.frame(freq(r_anthromes_yucatan))

df_freq <- merge(df_freq,df_anthromes_legend,by.x="value",by.y="Grid.Value")

df_freq$group_val <- c(1,1,2,2,2,3,3,3,3,3,4,4,4,5,5,6,6) #This reclasses anthromes categories in groups

r_group_anthromes_yucatan <- r_anthromes_yucatan
df_subs <- subset(df_freq,select=c("value","group_val"))
r_group_anthromes_yucatan <- subs(x=r_group_anthromes_yucatan,y=df_subs)

r_group_anthromes_yucatan_filename <- file.path(out_dir,paste0("r_group_anthromes_yucatan","_",out_suffix,".tif"))
writeRaster(r_group_anthromes_yucatan,r_group_anthromes_yucatan_filename,overwrite=T)

### Set color palette for the Anthromes group map
col_group <- c("dark grey","dark blue","yellowgreen","brown","darkgreen","darkorange")

plot(r_group_anthromes_yucatan,col=col_group,legend=FALSE,axes="FALSE")
legend("topleft",legend=unique(df_freq$GROUP),title="Anthromes",
       pt.cex=1.1,cex=1.2,fill=col_group,bty="n")

barplot(df_freq$count,names.arg=df_freq$LABEL,angle="90",las=2)
barplot(df_freq$count,names.arg=df_freq$group_val,angle="90",las=2)

######### Save combined figure 1 for paper: for fire and anthromes groups

#remove extreme values in fire count which is due to burning plant
minValue(r_fire_freq)
max_val <- maxValue(r_fire_freq)
r_fire_freq[r_fire_freq==max_val] <- NA

res_pix<-960
col_mfrow<-1
row_mfrow<-0.5 #set hight as half
png(filename=paste("Figure1_yucatan_fire_and_anthromes_at_10km2_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(1,2))

palette_fire <- colorRampPalette(c("lightgreen","yellow","red","darkred"))
#temp.colors <- matlab.like(18)
plot(r_fire_freq,main="Fire frequency count over 2001-2008",col=palette_fire(256) )

plot(r_group_anthromes_yucatan,col=col_group,legend=FALSE,axes="T",
     main="Anthromes Groups")
legend("topleft",legend=unique(df_freq$GROUP),title="Anthromes",
       pt.cex=1.1,cex=1.2,fill=col_group,bty="n")

#scale_position<-c(450000, 600000)
#arrow_position<-c(900000, 600000)
dev.off()

############## Model the relationship ########

r_c <-stack(r_anthromes_yucatan,r_group_anthromes_yucatan,r_fire_freq)
df_r_c <- as.data.frame(r_c)
names(df_r_c) <- c("cat","group_cat","fire_count")

df_r_c <- merge(df_r_c ,df_anthromes_legend,by.x="cat",by.y="Grid.Value")
dim(df_r_c) #2191x4

freq(r_anthromes_yucatan)

### Fix issue with Cropland groups and classes!!
df_r_c$GROUP <- as.character(df_r_c$GROUP)
table(df_r_c$GROUP)
#if(df_r_c)
unique(df_r_c$GROUP)

df_r_c$anth_cat <- as.factor(df_r_c$LABEL)
df_r_c$anth_group <-  as.factor(df_r_c$GROUP)

df_cat_mean <- aggregate(fire_count ~ cat, df_r_c, mean)
df_anthromes_legend
df_anth_cat_combined <- merge(df_cat_mean,df_anthromes_legend,by.x="cat",by.y="Grid.Value")
View(df_anth_cat_combined)

barplot(df_anth_cat_combined$fire_count,names.arg=df_anth_cat_combined$LABEL,angle="90",las=2)
barplot(df_anth_cat_combined$fire_count,names.arg=df_anth_cat_combined$LABEL,angle="90",las=2,horiz=T)

#boxplot(log_fire~cat, df_r_c, outline=F)
#,names=df_anthromes_legend$LABEL)

##############################
#### Try Poisson regression ###

df_r_c_all <- df_r_c
#df_r_c <- df_r_c_all
df_r_c <- df_r_c[!is.na(df_r_c$fire_count),]

sum(is.na(df_r_c$GROUP))

glm_fire <- glm(fire_count~cat,family="poisson",data=df_r_c)

glm_fire_anth_cat <- glm(fire_count~anth_cat,family="poisson",data=df_r_c)

glm_fire_anth_group <- glm(fire_count~anth_group,family="poisson",data=df_r_c)
summary(glm_fire_anth_group)

names(glm_nb_anth_group)
View(glm_nb_anth_cat$model)
df_r_c$fire_count
unique(df_r_c$anth_group)

#lm_fire_anth_cat <- lm(fire_count ~ anth_cat,df_r_c)
#summary(lm_fire_anth_cat)
#lm_fireby_group <- lm(fire_count ~ anth_cat,df_r_c)
#summary(lm_fire2)

histogram(log(df_r_c$fire_count))
df_r_c$log_fire <- log(df_r_c$fire_count)
hist(df_r_c$fire_count)

#> unique(df_r_c$cat)
#[1] NA 52 43 62 51 34 32 26 12 41 42 31 11 24 35 33 61 25
#> length(unique(df_r_c$cat))

glm_nb_fire_anth_group <- glm.nb(fire_count ~ anth_group, data = df_r_c)
glm_nb_fire_anth_cat <- glm.nb(fire_count ~ anth_cat, data = df_r_c)

agg_anth_group <- aggregate(fire_count ~anth_group,data = df_r_c,mean)
agg_anth_group$var <- (aggregate(fire_count ~anth_group,data = df_r_c,var))$fire_count
names(agg_anth_group) <- c("anth_group","mean","var")

View(agg_anth_group)

X2 <- 2 * (logLik(glm_nb_fire_anth_group) - logLik(glm_fire_anth_group))
X2
#'log Lik.' 135072.3 (df=8)
pchisq(X2, df = 1, lower.tail=FALSE)
#'log Lik.' 0 (df=8)

glm_nb_fire_anth_group$theta

#summary(anova(glm_nb_anth_group,test="LRT"))

comp_tukey_glm_fire_anth_group <- glht(glm_fire_anth_group,mcp(anth_group='Tukey'))

summary(comp_tukey_glm_fire_anth_group)

comp_tukey_glm_nb_fire_anth_group <- glht(glm_nb_fire_anth_group,mcp(anth_group='Tukey'))

## Generate summary Tukey output table for paper
#This is the ANOVA like analysis with pairwise comparison

summary(comp_tukey_glm_nb_fire_anth_group )
df_tukey_glm_nb_anth_group <- generate_tukey_table(comp_tukey_glm_nb_fire_anth_group ,format_digits,out_dir,out_suffix)

## write out table output:
df_tukey_glm_nb_anth_group_filename <- file.path(out_dir,paste0("df_tukey_glm_nb_anth_group_",out_suffix,".txt"))
write.table(df_tukey_glm_nb_anth_group, df_tukey_glm_nb_anth_group_filename,row.names=F,sep=",")

##### Now do model for rangeland and cropland classes

df_r_c_subset <- subset(df_r_c,anth_group%in%c("Croplands","Rangelands"))
df_r_c_subset$anth_cat <- as.character(df_r_c_subset$anth_cat)
df_r_c_subset$anth_cat <- as.factor(df_r_c_subset$anth_cat)
table(df_r_c_subset$anth_cat)

boxplot(fire_count ~ anth_cat,data=df_r_c_subset,outline=F,names,las=2)

glm_nb_fire_croplands_rangelands_anth_cat <- glm.nb(fire_count ~ anth_cat, data = df_r_c_subset)
summary(glm_nb_fire_croplands_rangelands_anth_cat)

comp_tukey_glm_nb_fire_croplands_rangelands_anth_cat<- glht(glm_nb_fire_croplands_rangelands_anth_cat,mcp(anth_cat='Tukey'))
summary(comp_tukey_glm_nb_fire_croplands_rangelands_anth_cat)

##### Now do model for cropland classes

df_r_c_subset <- subset(df_r_c,anth_group==("Croplands"))
df_r_c_subset$anth_cat <- as.character(df_r_c_subset$anth_cat)
df_r_c_subset$anth_cat <- as.factor(df_r_c_subset$anth_cat)
table(df_r_c_subset$anth_cat)

glm_nb_fire_croplands_anth_cat <- glm.nb(fire_count ~ anth_cat, data = df_r_c_subset)
summary(glm_nb_fire_croplands_anth_cat)

comp_tukey_glm_nb_fire_croplands_anth_cat<- glht(glm_nb_fire_croplands_anth_cat,mcp(anth_cat='Tukey'))
summary(comp_tukey_glm_nb_fire_croplands_anth_cat)

##### write out Tukey comparison
#### Output table for anth_cat for croplands

#debug(generate_tukey_table)
df_tukey_glm_nb_croplands_anth_cat <- generate_tukey_table(comp_tukey_glm_nb_fire_croplands_anth_cat,format_digits,out_dir,out_suffix)

## write out table output:
df_tukey_glm_nb_anth_group_filename <- file.path(out_dir,paste0("df_tukey_glm_nb_croplands_anth_cat_",out_suffix,".txt"))
write.table(df_tukey_glm_nb_croplands_anth_cat, df_tukey_glm_nb_anth_group_filename,row.names=F,sep=",")

############
#### Prepare figure 2 with boxplot ###

res_pix<-960
col_mfrow<-1
row_mfrow<-0.5 #set hight as half

png(filename=paste("Figure2_yucatan_fire_boxplots_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
par(mfrow=c(1,2))

#par(mar=c(6, 4.1, 4.1, 2.1))
par(mar=c(7, 4.5, 4.5, 2.3))

labels_group <- as.character(unique(df_r_c$anth_group))
labels_group <- c("Croplands ","Dense settlements ","Forested ","Rangelands ","Villages ","Wildlands ")
#labels_cat <- aas.character(unique(df_r_c_subset$anth_cat))
labels_cat <- c("Residential irrigated ","Residential rainfed mosaic ","Populated irrigated ",
                "Populated rainfed ","Remote croplands ")
#labels_cat <- c("Populated irrigated ","Populated rainfed ","Residential irrigated ","Residential rainfed mosaic ",,
#                ,"Remote cropland ")
#as.character(unique(df_r_c_subset$anth_cat))
boxplot(fire_count~anth_group, df_r_c, outline=F,main="Fire Frequency by Group Anthromes",
        col = "lightgray", xaxt = "n",  xlab = "")
# x axis with ticks but without labels
axis(1, labels = FALSE)

# Plot x labs at default x position
text(x =  seq_along(labels_group), y = par("usr")[3] - 1, srt = 40, adj = 1,
     labels = labels_group, xpd = TRUE)

boxplot(fire_count ~ anth_cat,data=df_r_c_subset,outline=F,main="Fire Frequency within Croplands Anthromes",
        col = "lightgray", xaxt = "n",  xlab = "")
# x axis with ticks but without labels
axis(1, labels = FALSE)

# Plot x labs at default x position
text(x =  seq_along(labels_cat), y = par("usr")[3] - 1, srt = 40, adj = 1,
     labels = labels_cat, xpd = TRUE)

#boxplot(count ~ spray, data = InsectSprays,
#        col = "lightgray", xaxt = "n",  xlab = "")

dev.off()

##################  END OF SCRIPT ######################

### relevant links:
#Poisson regression:
#https://onlinecourses.science.psu.edu/stat504/book/export/html/165
#http://www.ats.ucla.edu/stat/r/dae/poissonreg.htm
#https://onlinecourses.science.psu.edu/stat504/node/169
#Overdispersion test

#with(dat, tapply(daysabs, prog, function(x) {
#  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
#}))
#http://www.ats.ucla.edu/stat/r/dae/nbreg.htm
#Use the glm.nb from the MASS package to run the negative bionomial model
#### Check for over dispersion
#Aslo do Tukey comparison
#http://www.ats.ucla.edu/stat/r/dae/nbreg.htm
#http://stats.stackexchange.com/questions/71961/understanding-over-dispersion-as-it-relates-to-the-poisson-and-the-neg-binomial
#R> library(AER)
#R> data(RecreationDemand)
#R> rd <- glm(trips ~ ., data = RecreationDemand, family = poisson)
#R> dispersiontest(rd,trafo=1)

#### Check for over dispersion
#Aslo do Tukey comparison
#http://www.ats.ucla.edu/stat/r/dae/nbreg.htm
#http://stats.stackexchange.com/questions/71961/understanding-over-dispersion-as-it-relates-to-the-poisson-and-the-neg-binomial
#R> library(AER)
#R> data(RecreationDemand)
#R> rd <- glm(trips ~ ., data = RecreationDemand, family = poisson)
#R> dispersiontest(rd,trafo=1)

