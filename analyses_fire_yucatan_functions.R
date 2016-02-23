################################    FIRE YUCATAN PAPER  #######################################
############################      LAND COVER CHANGE          #######################################
#This script contains functions used for the analyses of Yucatan fire.
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#AUTHOR: Benoit Parmentier, Marco Millones                                                                      #
#DATE CREATED: 02/12/2016 
#DATE MODIFIED: 02/22/2016
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

run_multinom_mod <- function(list_models,model_type="multinom",y_var_name,data_df,ref_var_name=NULL){
  ##
  if(!is.null(ref_var_name)){
    data_df$y_var <- relevel(data_df[[y_var_name]], ref = ref_var_name)
  }else{
    data_df$y_var <- data_df[[y_var_name]]
  }
  list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!

  #formula_obj
  #run_multinom <- function(data_df,formula_obj){
  list_mod <- vector("list",length=length(list_formulas))
  if(model_type=="multinom"){
    #eval(list_formulas[[1]])
    #list_mod <- lapply(1:length(list_formulas),
    #                   FUN=function(i,list_formulas,data_df){formula_val <- list_formulas[[i]];multinom(formula_val,data=data_df)},
    #                  list_formulas=list_formulas,data_df=data_df)
    for(k in 1:length(list_formulas)){
      #formula_val <- list_formulas[[k]]
      formula_val <- list_models[[k]]
      mod <- multinom(as.formula(formula_val),data=data_df)#work around to solve issue in model object
      #mod <- multinom(y_var ~ cy_r, data=data_df)
      list_mod[[k]] <- mod
    }
    #mod <- multinom(formula_obj,data=data_df)
  }
  return(list_mod)
}
  
run_multinom_mod_mc <- function(i,list_param){
  #This function runs multinomial model using mclapply
  #
  #
  #Inputs:
  
  #Start script
  list_models <- list_param$list_models
  model_type <- list_param$model_type
  y_var_name <- list_param$y_var_name
  data_df <- list_param$data_df
  unique_val_y_var <- list_param$unique_val_y_var
  #ref_var_name <- list_param$ref_var_name
  out_suffix_s <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  ref_var_name <- unique_val_y_var[i]
  #out_suffix_s <- out_suffix_s[i]
  
  #test<- mclapply(1:length(unique_val_y_var),FUN=run_multinom_mod,list_models=list_models,
  #                model_type="multinom",y_var_name=y_var_name,data_df=data_df_spdf,ref_var_name=unique_val_y_var,
  #                mc.preschedule=FALSE,mc.cores = num_cores)
  
  #debug(run_multinom_mod)
  list_mod <- run_multinom_mod(list_models,model_type="multinom",y_var_name,data_df=data_df,ref_var_name=ref_var_name)
  names(list_mod) <- paste("ref_",ref_var_name,sep="")
  names_mod_obj <- file.path(".",paste("list_mod_","ref_",ref_var_name,"_",out_suffix_s,".RData",sep=""))
  save(list_mod,file= names_mod_obj)
  
  #for(i in 1:length(unique_val_y_var)){
  #  ref_var_name <- unique_val_y_var[i]
  #  #debug(run_multinom_mod)
  #  list_mod <- run_multinom_mod(list_models,model_type="multinom",y_var_name,data_df=data_df_spdf,ref_var_name=ref_var_name)
  #   names(list_mod) <- paste("ref_",ref_var_name)
  # names_mod_obj <- file.path(".",paste("list_mod_","ref_",ref_var_name,"_",out_suffix_s,".RData",sep=""))
  # save(list_mod,file= names_mod_obj)
  
  #  list_mod_obj[[i]]<- list_mod
  #}
  #
  return(list_mod)
}

multinomial_model_fun<-function(list_models,model_type="multinom",y_var_name,data_df=data_df_spdf,ref_var_name,zonal_var_name,num_cores,out_suffix_s,out_dir){
  #This function runs multinomial models 
  #Inputs:
  #1)list_models: list of models formul as string/character
  #2)model_type: currently only "multinom"
  #3)y_var_name: name of the variable to use as dependent variable
  #4)data_df: data.frame containing the variables
  #5)ref_var_name: reference name for the multinomial model
  #6)zonal_var_name: stratum used to run models as subset
  #7)num_cores: number of cores to use
  #8)out_suffix_s: output suffix
  #9)out_dir: output dir used in the work
  #Output:
  #
  
  ### Start of script ###
  
  #unique_val_zones <- (unique(data_df_spdf[[zonal_var_name]]))
  #unique_val_zones <- unique_val_zones[!is.na(unique_val_zones)]
  unique_val_y_var <- as.numeric(unique(data_df[[y_var_name]]))
  unique_val_y_var <- sort(unique_val_y_var[!is.na(unique_val_y_var)]) #values 1,2,3
  
  #
  #
  
  #list_mod_obj <- vector("list",length=length(unique_val_y_var))
  #list_models,model_type="multinom",y_var_name,data_df,ref_var_name=NULL
  
  #out_suffix_s <- 
  list_param_multinom <- list(list_models,model_type,y_var_name,data_df,unique_val_y_var,out_suffix_s,out_dir)
  names(list_param_multinom) <- c("list_models","model_type","y_var_name","data_df","unique_val_y_var","out_suffix","out_dir")
  
  #debug(run_multinom_mod_mc)
  #test<- run_multinom_mod_mc(1,list_param_multinom)
  model_obj <- mclapply(1:length(unique_val_y_var),FUN=run_multinom_mod_mc,list_param=list_param_multinom,
                        mc.preschedule=FALSE,mc.cores = num_cores)
  #r_var_s <- mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  
  return(model_obj)
} 

extract_coef_p_values <- function(mod){
  
  summary_mod <- summary(mod)
  z <- summary_mod$coefficients/summary_mod$standard.errors
  #2-tailed z test
  p <- (1 - pnorm(abs(z), 0, 1))*2 #95% interval Wald
  
  summary_coefficients <- summary_mod$coefficients
  summary_standard_errors <- summary_mod$standard.errors
  
  summary_obj <- list(p,z,summary_coefficients,summary_standard_errors)
  names(summary_obj) <- c("p","z","summary_coefficients","summary_standard_errors")
  return(summary_obj)
}

extract_multinom_mod_information <- function(mod){
  #modify to use mclapply later!!!
  
  
  ##### Start script ###
  
  if(class(mod)!="list"){
    list_mod <- list(mod)
  }else{
    list_mod <- mod
    rm(mod)
  }
  
  names_mod <- paste("mod",1:length(list_mod),sep="")
  AIC_values <- unlist(lapply(list_mod,function(x){x$AIC}))
  names(AIC_values) <- names_mod
  list_coef <- lapply(list_mod,function(x){summary(x)$coefficients})
  #list_formulas <- lapply(list_mod,function(x){summary(x)$formula})
  list_extract_coef_p_values <- lapply(list_mod,FUN=extract_coef_p_values)
  names(list_extract_coef_p_values) <-names_mod
  multinom_extract_obj <- list(AIC_values,list_coef,list_extract_coef_p_values)
  names(multinom_extract_obj) <- c("AIC_values","list_coef","list_extract_coef_p_values")
  return(multinom_extract_obj)
}

### Add function to reformat the tables here...


############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
#http://data.princeton.edu/wws509/r/c6s2.html

