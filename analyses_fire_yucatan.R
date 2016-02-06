################################    FIRE YUCATAN PAPER  #######################################
############################      LAND COVER CHANGE          #######################################
#This script reads data extracted from MODIS FIRE time series                             
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#AUTHOR: Benoit Parmentier, Marco Millones                                                                      #
#DATE CREATED: 02/06/2016 
#DATE MODIFIED: 02/06/2016
#Version: 1
#PROJECT: Land cover Change Yucatan, Marco Millones
#   
#COMMENTS: Initial commit
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

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_05172015_functions.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

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
data_CI_fname <- "/home/bparmentier/Google Drive/FireYuca_2016/yucatan_fire.xlsx" #contains the whole dataset
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
data_CI <-read.table(data_CI_fname, sheet=1)

#write.table(data_Hansen,file.path(in_dir,"data_Hansen_02062016.txt"),sep=",")
write.table(data_CI,file.path(in_dir,"data_CI_02062016.txt"),sep=",")

data_CI_spdf <- data_CI
coordinates(data_CI_spdf) <- coord_names

### Quick check of plotting variables
spplot(data_CI_spdf,"cattledensity")
spplot(data_CI_spdf,"FIRE_pre07")
spplot(data_CI_spdf,"dist_roads")

#Create the count variable
fire_modis_col <- c("FIRE_2000","FIRE_2001","FIRE_2002","FIRE_2003","FIRE_2004","FIRE_2005","FIRE_2006","FIRE_2007")
data_CI$FIRE_freq <- colSums(data_CI[,fire_modis_col])
data_CI$FIRE_intensity <- data_CI$fire_freq/8
data_CI$FIRE_bool <- data_CI$FIRE_freq > 0

#### Run modeling series A ###

#run overall but also state by state
mod1a <- multinom(fpnfpch ~ cy_q , data = data_CI)
mod2a <- multinom(fpnfpch ~ cy_q + FIRE_pre07, data = data_CI)
mod3a <- multinom(fpnfpch ~ cy_q + FIRE_intensity, data = data_CI)
mod4a <- multinom(fpnfpch ~ cy_q + FIRE_bool, data = data_CI)

mod5a <- multinom(fpnfpch ~ cy_q +  
                            dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                            elevation + slope_deg + landcover + cattledensity + ejido + 
                            popdens_change + precip + protegidas + soil
                            , data = data_CI)

mod6a <- multinom(fpnfpch ~ cy_q + FIRE_pre07 + 
                    dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                    elevation + slope_deg + landcover + cattledensity + ejido + 
                    popdens_change + precip + protegidas + soil
                  , data = data_CI)

mod7a <- multinom(fpnfpch ~ cy_q + FIRE_intensity + 
                    dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                    elevation + slope_deg + landcover + cattledensity + ejido + 
                    popdens_change + precip + protegidas + soil
                  , data = data_CI)

mod8a <- multinom(fpnfpch ~ cy_q + FIRE_freq + 
                    dist_disturbed + dist_loc + dist_rur + dist_urb + dist_roads + 
                    elevation + slope_deg + landcover + cattledensity + ejido + 
                    popdens_change + precip + protegidas + soil
                  , data = data_CI)

#extract the AIC and formula!!!
#extract slope of fire?
#
list_moda <- list(mod1a,mod2a,mod3a,mod4a,mod5a,mod6a,mod7a,mod8a)

############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
