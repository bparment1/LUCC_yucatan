#################################    FIRE YUCATAN PAPER  #######################################
#############################         LAND COVER CHANGE          #######################################
#This script reads data extracted from MODIS FIRE time series                             
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#
#AUTHOR: Benoit Parmentier                                                                #
#DATE CREATED: 02/06/2016 
#DATE MODIFIED: 03/09/2017
#Version: 1
#PROJECT: Land cover Change Yucatan with Marco Millones 
#   
#COMMENTS:  
#TODO:
#
#COMMITS: rerunning models for study/paper
#
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

###### Functions used in this script

functions_analyses_script <- "analyses_fire_yucatan_functions_05282016.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/FireYuca_2016/R_scripts" #path to script #PARAM 2
source(file.path(script_path,functions_analyses_script)) #source all functions used in this script 1.

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
out_suffix <-"yucatan_CI_analyses_03092017" #output suffix for the files and ouptu folder #PARAM 8
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
run_relevel <- TRUE
num_cores <- 3

list_models<-c("y_var ~ cy_r",
               "y_var ~ cy_r + FIRE_pre07",
               "y_var ~ cy_r + FIRE_intensity",
               "y_var ~ cy_r + FIRE_bool",
               "y_var ~ cy_r +  
                      dist_disturbed + dist_NF00 + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_r + FIRE_bool +
                      dist_disturbed + dist_NF00 + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_r + FIRE_pre07 + 
                      dist_disturbed + dist_NF00 + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_r + FIRE_intensity + 
                      dist_disturbed + dist_NF00 + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil",
               "y_var ~ cy_r + FIRE_freq + 
                      dist_disturbed + dist_NF00 + dist_rur + dist_urb + dist_roads + 
                      elevation + slope_deg + cattledensity + ejido + 
                      popdens_change + precip + protegidas + soil"
)

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

data_df <-read.table(data_CI_fname,stringsAsFactors=F,header=T,sep=",")
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
#> sum(table(data_df_spdf$state))
#[1] 136311
#> nrow(data_df_spdf)
#[1] 136904

### Quick check of plotting variables
spplot(data_df_spdf,"cattledensity")
spplot(data_df_spdf,"FIRE_pre07")
spplot(data_df_spdf,"dist_roads")
spplot(data_df_spdf,"state")

hist(data_df[[y_var_name]])
table(data_df[[y_var_name]])
barplot(table(data_df[[y_var_name]]))

#### Part II: run modeling ####

#### Run modeling series A ###

#run overall for the whole Yucatan peninsula

#make this a function?

#debug(run_multinom_mod)
if(run_relevel==TRUE){
  data_df_spdf[[y_var_name]] <- as.factor(data_df_spdf[[y_var_name]])
}

#Run for the overall state
out_suffix_s <- paste("overall_",out_suffix,sep="")
#debug(multinomial_model_fun)

data_df_overall_model_obj <- multinomial_model_fun(list_models,model_type="multinom",y_var_name,data_df=data_df_spdf,ref_var_name,zonal_var_name,num_cores,out_suffix_s,out_dir)
names_mod_obj <- file.path(out_dir,paste("data_df_overall_model_obj_",out_suffix_s,".RData",sep=""))
save(data_df_overall_model_obj,file= names_mod_obj)

#### Run modeling series A ###

#run state by state
unique_val_zones <- (unique(data_df_spdf[[zonal_var_name]]))
unique_val_zones <- unique_val_zones[!is.na(unique_val_zones)]

#> unique_val_zones
#[1] "Yucatan"      "Quintana Roo" "Campeche"
list_by_region_model_obj <- vector("list",length=length(unique_val_zones))
for(i in 1:length(unique_val_zones)){
 
  val_zone <- paste(unlist(strsplit(unique_val_zones[i]," ")),collapse = "_") #e.g. state here
  out_suffix_s <- paste(val_zone,out_suffix,sep="_") #use data by state
  data_subset <- subset(data_df_spdf,data_df_spdf$state==unique_val_zones[i])
  data_df_by_region_model_obj <- multinomial_model_fun(list_models,model_type="multinom",y_var_name,data_df=data_subset,ref_var_name,zonal_var_name,num_cores,out_suffix_s,out_dir)
  
  #names(list_mod) <- paste("val_zone_",unique_val_zones[i],sep="")
  names_mod_obj <- file.path(out_dir,paste("data_df_by_region_model_obj_","val_zone_","_",out_suffix_s,".RData",sep=""))
  save(data_df_by_region_model_obj,file= names_mod_obj)
  list_by_region_model_obj[[i]] <- data_df_by_region_model_obj
  names(list_by_region_model_obj[[i]]) <- val_zone
}

#### PART III: extract information from models ######


list_model_objects <- list.files(out_dir,"*model_obj*",full.names=T)
#debug(extract_multinom_mod_information)  
ref_var <- c(1,2,3)
list_region_name <- c("Campeche","Quintana_Roo","Yucatan","overall")
region_name <- list_region_name

### Make this part a function:
#debug(generate_summary_tables_from_models)
test <- generate_summary_tables_from_models(list_model_objects,ref_var,region_name,out_suffix,out_dir)

### Get accuracy information for models: AIC, deviance, log likelihood, nb of observation etc.
#debug(extract_model_metrics_accuracy)
                                            
tb_models_accuracy <- extract_model_metrics_accuracy(list_model_objects,region_name,num_cores,out_suffix,out_dir)
#list_extract_mod_yucatan_CI_analyses_04192016.RData

## Get odds ratio

#### Add extraction of fitted values and residuals??

#length(mod$fitted.values)
#dim(mod$fitted.values)

##### PART IV: Get general information about the data ######

xtabs_tb <- xtabs(~fpnfpch+FIRE_bool,data_df)
attach(data_df)
mytable <- table(fpnfpch,FIRE_bool) # A will be rows, B will be columns 
mytable # print table 

crosstab(data_df, row.vars = "FIRE_bool", col.vars = "fpnfpch", type = "f")
CrossTable(data_df$fpnfpch, data_df$FIRE_bool)

data_df$FIRE_bool <- as.numeric(data_df$FIRE_bool)


############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
#http://data.princeton.edu/wws509/r/c6s2.html

#> mod2a$coefnames
#[1] "(Intercept)" "cy_q"        "FIRE_pre07" 
#> coefficients(mod2a)
#(Intercept)      cy_q FIRE_pre07
#2  -0.1389069 -0.382720  0.2230844
#3  -2.2844927 -3.746468  0.1316694

#z <- summary(mod8a)$coefficients/summary(mod8a)$standard.errors
#z

#2-tailed z test
#p <- (1 - pnorm(abs(z), 0, 1))*2
#p



#
# Browse[4]> extract_mod[[i]]
# $p
# (Intercept) cy_r
# 1           0    0
# 2           0    0
# 
# $z
# (Intercept)     cy_r
# 1    9.603259 16.92589
# 2   13.327639 12.93944
# 
# $summary_coefficients
# (Intercept)     cy_r
# 1    1.835426 5.676169
# 2    2.535650 4.330786
# 
# $summary_standard_errors
# (Intercept)      cy_r
# 1   0.1911253 0.3353542
# 2   0.1902550 0.3346966

