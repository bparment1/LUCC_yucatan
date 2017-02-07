#################################    FIRE YUCATAN PAPER  #######################################
#############################         LAND COVER CHANGE          #######################################
#This a function script for Yucatan fire modeling.
#The script reads data extracted from MODIS FIRE time series                             
#The goal is to show how fire can be used as proxy for land cover change in the region.                                                     #
#The link between fire as a tool for anthropogenic changes is also explore with the anthrome dataset.
#We used a non negative binomial model with Tukey to examine differences between anthromes classes.
#
#AUTHOR: Benoit Parmentier                                                                #
#DATE CREATED: 02/06/2016 
#DATE MODIFIED: 02/07/2017
#Version: 1
#PROJECT: Land cover Change Yucatan with Marco Millones 
#   
#COMMIT: modifying Tukey table and some clean up of code
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
library(parallel)                        # parralel programming and multi cores
library(plyr)
library(rgeos)                           # topology and vector spatial queries and operations
library(afex)                            # functions related to ANOVA
library(car)
library(MASS)                            # contains negative binomial model
library(multcomp)                        # contains Tukey comparison

###### Functions used in this script

generate_tukey_table <- function(comp_tukey_glm,format_digits,out_dir,out_suffix){
  #Process output from glht Tukey pairwise comparison
  #This assumes Poisson or Negative binomial model
  #
  ## Inputs:
  #comp_tukey_glm
  #
  # Outputs
  #
  
  ###########################
  
  #### Begin Script ###
  
  summary_comp_tukey <- summary(comp_tukey_glm)
  pairwise_comparison <- names(summary_comp_tukey$test$coefficients)
  df_tukey <- as.data.frame(
    cbind(pairwise_comparison=pairwise_comparison,
          format(summary_comp_tukey$test$coefficients,digits=format_digits),
          format(summary_comp_tukey$test$sigma,digits=format_digits),
          format(summary_comp_tukey$test$tstat,digits=format_digits),
          format(as.numeric(summary_comp_tukey$test$pvalues,digits=format_digits)))
  )
  
  names(df_tukey) <- c("pairwise_comparison","coefficients","standard_error","t-stat","p_values")
  row.names(df_tukey) <- NULL
  df_tukey$p_values <- round(as.numeric(as.character(df_tukey$p_values)), digits = 4)
  
  return(df_tukey)
  
}

############### END OF SCRIPT ###################

