################################    FIRE YUCATAN PAPER  #######################################
############################      LAND COVER CHANGE          #######################################
#This script contains functions used for the analyses of Yucatan fire.
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#AUTHOR: Benoit Parmentier, Marco Millones                                                                      #
#DATE CREATED: 02/12/2016 
#DATE MODIFIED: 02/12/2016
#Version: 1
#PROJECT: Land cover Change Yucatan, Marco Millones
#   
#COMMENTS: Initial commit: Separation between function script and main script
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

run_multinom_mod <- function(list_models,model_type="multinom",y_var_name,data_df){
  ##
  
  list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!
  data_df$y_var <- data_df[[y_var_name]]
  #formula_obj
  #run_multinom <- function(data_df,formula_obj){
  if(model_type=="multinom"){
    list_mod <- lapply(1:length(list_formulas),
                       FUN=function(i,list_formulas,data_df){multinom(list_formulas[[i]],data=data_df)},
                       list_formulas=list_formulas,data_df=data_df)
    #mod <- multinom(formula_obj,data=data_df)
  }
  return(list_mod)
}
  
extraction_of_information <- function(list_mod){
  
  #list_moda <- list(mod1a,mod2a,mod3a,mod4a,mod5a,mod6a,mod7a,mod8a)
  AIC_values <- unlist(lapply(list_mod,function(x){x$AIC}))
  list_coef <- lapply(list_mod,function(x){summary(x)$coefficients})
  list_formulas <- lapply(list_mod,function(x){summary(x)$formula})
  multinom_extract_obj <- list(AIC_values,list_coef,list_formulas)
  names(multinom_extract_obj) <- c("AIC_values","list_coef","list_formulas")
  return(multinom_extract_obj)
}

############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
#http://data.princeton.edu/wws509/r/c6s2.html

