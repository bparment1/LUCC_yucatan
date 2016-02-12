################################    FIRE YUCATAN PAPER  #######################################
############################      LAND COVER CHANGE          #######################################
#This script reads data extracted from MODIS FIRE time series                             
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#AUTHOR: Benoit Parmentier, Marco Millones                                                                      #
#DATE CREATED: 02/06/2016 
#DATE MODIFIED: 02/12/2016
#Version: 1
#PROJECT: Land cover Change Yucatan, Marco Millones
#   
#COMMENTS:  Separation between function script and main script
#TODO:

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

###### Functions used in this script

functions_analyses_script <- "analyses_fire_yucatan_functions_02122016.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/FireYuca_2016/R_scripts" #path to script #PARAM 2
source(file.path(script_path,functions_analyses_script)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir<-"/home/bparmentier/Google Drive/FireYuca_2016/"

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"yucatan_analyses_02062016" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

data_Hansen_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/Hansen_fire.xlsx" #contains the whole dataset
#data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/yucatan_fire.xlsx" #contains the whole dataset
data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/yucatan_fire.csv" #contains the whole dataset
#data_CI_02062016.txt
coord_names <- c("POINT_X","POINT_Y") #PARAM 11

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#data_tb <-read.table(data_fname,sep=",",header=T)

#Import data from excel
#data_Hansen <-read.xls(data_Hansen_fname, sheet=1)
#data_CI <-read.xls(data_CI_fname, sheet=1)
#data_Hansen <-read.xls(data_Hansen_fname, sheet=1)
data_CI <-read.table(data_CI_fname,stringsAsFactors=F,header=T,sep=",")

#write.table(data_Hansen,file.path(in_dir,"data_Hansen_02062016.txt"),sep=",")
#write.table(data_CI,file.path(in_dir,"data_CI_02062016.txt"),header=T,sep=",")

data_CI_spdf <- data_CI
coordinates(data_CI_spdf) <- coord_names

### Quick check of plotting variables
spplot(data_CI_spdf,"cattledensity")
spplot(data_CI_spdf,"FIRE_pre07")
spplot(data_CI_spdf,"dist_roads")

#Create the count variable
fire_modis_col <- c("FIRE_2000","FIRE_2001","FIRE_2002","FIRE_2003","FIRE_2004","FIRE_2005","FIRE_2006","FIRE_2007")
data_CI$FIRE_freq <- colSums(data_CI[,fire_modis_col])
data_CI$FIRE_intensity <- data_CI$FIRE_freq/8
data_CI$FIRE_bool <- data_CI$FIRE_freq > 0

hist(data_CI$fpnfpch)
table(data_CI$fpnfpch)

#### Run modeling series A ###

#run overall but also state by state

#make this a function?

list_models<-c("y_var ~ cy_q",
               "y_var ~ cy_q + FIRE_pre07",
               "y_var ~ cy_q + FIRE_intensity",
               "y_var ~ cy_q + FIRE_bool",
               "y_var ~ cy_q +  
                      dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + landcover + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_q + FIRE_pre07 + 
                      dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + landcover + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_q + FIRE_intensity + 
                      dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + landcover + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_q + FIRE_freq + 
                      dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + landcover + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil"
               )

y_var_name <- "fpnfpch"
debug(run_multinom_mod)
list_mod <- run_multinom_mod(list_models,model_type="multinom",y_var_name ,data_df=data_CI)

  
# mod1a <- multinom(fpnfpch ~ cy_q , data = data_CI)
# mod2a <- multinom(fpnfpch ~ cy_q + FIRE_pre07, data = data_CI)
# mod3a <- multinom(fpnfpch ~ cy_q + FIRE_intensity, data = data_CI)
# mod4a <- multinom(fpnfpch ~ cy_q + FIRE_bool, data = data_CI)
#   
# mod5a <- multinom(fpnfpch ~ cy_q +  
#                       dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
#                       elevation + slope_deg + landcover + cattledensity + ejido + 
#                       popdens_change + precip + protegidas + soil
#                     , data = data_CI)
#   
# mod6a <- multinom(fpnfpch ~ cy_q + FIRE_pre07 + 
#                       dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
#                       elevation + slope_deg + landcover + cattledensity + ejido + 
#                       popdens_change + precip + protegidas + soil
#                     , data = data_CI)
#   
# mod7a <- multinom(fpnfpch ~ cy_q + FIRE_intensity + 
#                       dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
#                       elevation + slope_deg + landcover + cattledensity + ejido + 
#                       popdens_change + precip + protegidas + soil
#                     , data = data_CI)
#   
# mod8a <- multinom(fpnfpch ~ cy_q + FIRE_freq + 
#                       dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
#                       elevation + slope_deg + landcover + cattledensity + ejido + 
#                       popdens_change + precip + protegidas + soil
#                     , data = data_CI)
  
mod_names <-  c("mod1a","mod2a","mod3a","mod4a","mod5a","mod6a","mod7a","mod8a")
df_mod <- data.frame(mod_names=mod_names,AIC=AIC_values)
barplot(df_mod$AIC,names="mod_names",names.arg=mod_names)

##Note that multinom does not calculate p-values but we can carry out a Wald test...

#> mod2a$coefnames
#[1] "(Intercept)" "cy_q"        "FIRE_pre07" 
#> coefficients(mod2a)
#(Intercept)      cy_q FIRE_pre07
#2  -0.1389069 -0.382720  0.2230844
#3  -2.2844927 -3.746468  0.1316694

z <- summary(mod8a)$coefficients/summary(mod8a)$standard.errors
z

#2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1))*2
p

############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
#http://data.princeton.edu/wws509/r/c6s2.html

