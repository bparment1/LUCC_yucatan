################################    FIRE YUCATAN PAPER  #######################################
############################      LAND COVER CHANGE          #######################################
#This script contains functions used for the analyses of Yucatan fire.
#The goal is to show how fire can be used as proxy for land cover change in the region.      
#This script implements a spatial logit model with explantory variables.                                                     #
#
#AUTHOR: Benoit Parmentier                                                                    #
#DATE CREATED: 02/12/2016 
#DATE MODIFIED: 05/28/2016
#Version: 1
#PROJECT: Land cover Change Yucatan with Marco Millones
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
  ## extract the coefficients from the model and exponentiate
  #exp(coef(test))  
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


## extract the coefficients from the model and exponentiate
#exp(coef(test))


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
  names(list_coef) <- names_mod
  #adding extraction of odds ratio
  list_odds <- lapply(list_mod,function(x){exp(coef(x))})
  names(list_odds) <- names_mod
  #list_formulas <- lapply(list_mod,function(x){summary(x)$formula})
  list_extract_coef_p_values <- lapply(list_mod,FUN=extract_coef_p_values)
  names(list_extract_coef_p_values) <- names_mod
  multinom_extract_obj <- list(AIC_values,list_coef,list_extract_coef_p_values,list_odds)
  names(multinom_extract_obj) <- c("AIC_values","list_coef","list_extract_coef_p_values","list_odds")
  return(multinom_extract_obj)
}

### Add function to reformat the tables here...

create_summary_tables_reg_multinom <- function(list_extract_mod,reg_param,ref_var,list_region_name,list_models){
  #
  #Inputs:
  #list_extract_mod
  #reg_param,
  #ref_var
  #list_region_name: the name for each region e.g. Campeche, Yucatan,Quitano Ro, Overall
  #list_models
  
  list_df_table_mod <- vector("list",length(list_extract_mod))
  ## if list_extract_mod is a list greater than 1 then exctract for each.
  ## in this context, the length is 4: Campeche, Yucatan,Quitano Ro, Overall
  for(i in 1:length(list_extract_mod)){
    #reg_param <- list_tables_reg_param[1] #start with p value first
    #list_models
    #region_name <- names(list_extract_mod)[i]
    region_name <- list_region_name[i]
    list_df_table_ref <- vector("list",length(ref_var))
    #extract_mod <- list_extract_mod[[4]][[1]]$list_extract_coef_p_values
    for(k in 1:length(ref_var)){
      extract_mod <- list_extract_mod[[i]][[k]]$list_extract_coef_p_values #this is the general table/list
      ref_var_tmp <- k
      #debug(produce_regression_parameter_table)
      test_df_tmp <- produce_regression_parameter_table(reg_param,extract_mod,ref_var_tmp,list_models,region_name)
      list_df_table_ref[[k]] <- test_df_tmp
    }
    names(list_df_table_ref) <- as.character(ref_var_tmp)
    list_df_table_mod[[i]]<-list_df_table_ref
  }
  names(list_df_table_mod) <- list_region_name
  return(list_df_table_mod)
}

produce_regression_parameter_table <- function(reg_param,extract_mod,ref_var_tmp,list_models,region_name){
  #This function reorganizes parameters extracted from multinom.
  #Parameters are, p, z, coefficients and standard errors and odds
  #Inputs:
  #1)reg_param <- regression parameters, p, z, coeff, standard errors, odds 
  #2)extract_mod <- extact object from function extract
  #3)ref_var_tmp: value for the reference variable, in this context 1,2,3
  #4)list_models: model formulas as string/char
  #5)region_name: overall, Campeche, Yucatan, Quintana Roo
  #
  
  list_df_tmp <- vector("list",length=length(list_models))
  #loop through models...
  for(i in 1:length(list_models)){
    if(reg_param=="odds"){
      reg_param_tmp <- "summary_coefficients"
      df_tmp <- as.data.frame((extract_mod[[i]][[reg_param_tmp]])) #e.g. p variable
      df_tmp <- exp(df_tmp)
    }else{
      df_tmp <- as.data.frame((extract_mod[[i]][[reg_param]])) #e.g. p variable
    }

    df_tmp$ref_eqt <- rownames(df_tmp)
    df_tmp$ref_var <- ref_var_tmp
    df_tmp$mod_name <- paste("mod",i,sep="")
    df_tmp$region <- region_name
    list_df_tmp[[i]] <- df_tmp
  }
  
  test_df <- do.call(rbind.fill,list_df_tmp)
  
  return(test_df)
}

generate_summary_tables_from_models <- function(list_models_objects,ref_var,region_name,out_suffix_s,out_dir){
  #This function generate summary tables from multinom model produced beforehand.
  #When a list of models is provided csv files are produced with each models term for p, coefficients,
  #standard errors and z values of each term.
  #
  #Inputs
  #1) list_models_objects
  #2) ref_var
  #3) region_name
  #4) out_suffix_s
  #5) out_dir
  #Outputs
  #
  
  ## Functions used in the script ####
  
  write_table_fun <- function(i,list_tb_summary,out_filename){
    #write.table(tb_summary,file=paste("tb_summary_p_Campeche_",out_suffix,".txt",sep=""),sep=",")
    write.table(list_tb_summary[[i]],file=out_filename[i],sep=",")
    return(out_filename[i])
  }
  
  ## Start ##
  
  list_extract_mod <- vector("list",length=length(list_model_objects))
  names(list_extract_mod) <- region_name 
  
  for(i in 1:length(list_model_objects)){
    #Reading object files produced earlier:
    model_obj <- load_obj(list_model_objects[[i]]) #ref 1 for overall model?
    #undebug(extract_multinom_mod_information)
    #test2 <- extract_multinom_mod_information(model_obj[[1]])
    #
    list_extract_mod[[i]] <- lapply(model_obj,FUN=extract_multinom_mod_information)
  }
  list_extract_mod_fname <- paste("list_extract_mod_",out_suffix_s,".RData",sep="")
  save(list_extract_mod,file=list_extract_mod_fname)
  #test <- extract_multinom_mod_information(model_obj[[1]])#ref 1
  
  names(list_extract_mod[[1]][[1]]) #
   
  #names(list_extract_mod[[1]][[1]]) #
  #[1] "AIC_values"                 "list_coef"                  "list_extract_coef_p_values"
  #[4] "list_odds" 
  
  #list_extract_mod <- list_extrat_mod
  
  #list_extract_mod[[4]][[1]]$list_extract_coef_p_values$mod1$p #this is for overall
  #list_extract_mod[[4]][[1]]$list_extract_coef_p_values$mod2$p
  #list_tables_reg_param <- names(list_extract_mod[[4]][[1]]$list_extract_coef_p_values$mod2)
  
  #> names(list_extract_mod[[4]][[1]]$list_extract_coef_p_values$mod2)
  #[1] "p"                       "z"                       "summary_coefficients"    "summary_standard_errors"
  
  #list_region_name <- c("Campeche","Quintana_Roo","Yucatan","overall")
  
  ### Generate summary tables with p, z, coef, std_errors for each term and model by region
  
  #undebug(create_summary_tables_reg_multinom)
  reg_param<- "p"
  tb_summary_p <- create_summary_tables_reg_multinom(list_extract_mod,reg_param,ref_var,list_region_name,list_models)
  reg_param<- "summary_coefficients"
  tb_summary_coef <- create_summary_tables_reg_multinom(list_extract_mod,reg_param,ref_var,list_region_name,list_models)
  reg_param<- "summary_standard_errors"
  tb_summary_std_errors <- create_summary_tables_reg_multinom(list_extract_mod,reg_param,ref_var,list_region_name,list_models)
  reg_param<- "z"
  tb_summary_z <- create_summary_tables_reg_multinom(list_extract_mod,reg_param,ref_var,list_region_name,list_models)
  
  ###Now extract odds
  #reg_param <- "odds"
  #(list_extract_mod[[1]][[2]])$list_odds
  #lapply(list_extract_mod,Fun=)
  #(list_extract_mod[[1]][[1]])$list_odds
  reg_param <- "odds"
  tb_summary_odds <- create_summary_tables_reg_multinom(list_extract_mod,reg_param,ref_var,list_region_name,list_models)
  
  ##Collapse tables for each region
  
  list_tb_summary_p <- lapply(tb_summary_p, function(x){do.call(rbind,x)}) #equal to region's length
  list_tb_summary_coef <- lapply(tb_summary_coef, function(x){do.call(rbind,x)})
  list_tb_summary_std_errors <- lapply(tb_summary_std_errors, function(x){do.call(rbind,x)})
  list_tb_summary_z <- lapply(tb_summary_z, function(x){do.call(rbind,x)})
  list_tb_summary_odds <- lapply(tb_summary_odds, function(x){do.call(rbind,x)})
  
  list_out_suffix_s <- paste(region_name,out_suffix_s,sep="_")
  
  out_filename <- file.path(out_dir,paste("tb_summary_p_",list_out_suffix_s,".txt",sep=""))
  out_filename_summary_p <- lapply(1:length(list_tb_summary_p),FUN=write_table_fun,list_tb_summary=list_tb_summary_p,out_filename=out_filename)
  out_filename <- file.path(out_dir,paste("tb_summary_coef_",list_out_suffix_s,".txt",sep=""))
  out_filename_summary_coef <- lapply(1:length(list_tb_summary_coef),FUN=write_table_fun,list_tb_summary=list_tb_summary_coef,out_filename=out_filename)
  out_filename <- file.path(out_dir,paste("tb_summary_std_errors_",list_out_suffix_s,".txt",sep=""))
  out_filename_summary_std_errors <- lapply(1:length(list_tb_summary_p),FUN=write_table_fun,list_tb_summary=list_tb_summary_std_errors,out_filename=out_filename)
  out_filename <- file.path(out_dir,paste("tb_summary_z_",list_out_suffix_s,".txt",sep=""))
  out_filename_summary_z <- lapply(1:length(list_tb_summary_z),FUN=write_table_fun,list_tb_summary=list_tb_summary_z,out_filename=out_filename)
  
  out_filename <- file.path(out_dir,paste("tb_summary_odds_",list_out_suffix_s,".txt",sep=""))
  out_filename_summary_odds <- lapply(1:length(list_tb_summary_odds),FUN=write_table_fun,list_tb_summary=list_tb_summary_odds,out_filename=out_filename)
  
  ###### Prepare output object

  list_out_filename <- list(out_filename_summary_p,out_filename_summary_coef,out_filename_summary_std_errors,
                            out_filename_summary_z,out_filename_summary_odds)
  list_out_filename <- unlist(list_out_filename)
  
  ##Now gather AIC and 
  
    
  ### Prepare objects to return
  
  
  return(list_out_filename)
}


### Now get AIC and log likelihood
#extract_model_metrics_accuracy <- function(list_extract_mod,out_suffix,out_dir){
extract_model_metrics_accuracy <- function(list_model_objects,region_name,num_cores,out_suffix,out_dir){
  #Extract AIC, loglikelihood
  #Extract odd ratio and interval
  
  ##use parallelization...
  list_df_val <- vector("list",length=length(list_model_objects))
  
  for(i in 1:length(list_model_objects)){
    #Reading object files produced earlier:
    model_obj <- load_obj(list_model_objects[[i]]) #ref 1 for overall model?
    #undebug(extract_multinom_mod_information)
    #list_extract_mod[[i]] <- lapply(model_obj,FUN=extract_multinom_mod_information)
    #list_extract_mod[[i]]
    list_mod <- model_obj[[1]] #use ref1
    #mod <- list_mod[[1]]
    names_mod <- paste("mod",1:length(list_mod),sep="")
    AIC_values <- unlist(lapply(list_mod,function(x){x$AIC}))
    res_deviance_values <- unlist(lapply(list_mod,function(x){x$deviance}))
    names(res_deviance_values) <- names_mod
    loglikelihood_values <- unlist(lapply(list_mod,function(x){x$value})) #this value is about 1/2 of deviance
    names(loglikelihood_values) <- names_mod
    n_obs <- unlist(lapply(list_mod,function(x){nrow(x$fitted.values)})) #this value is about 1/2 of deviance
    names(n_obs) <- names_mod    
    
    ##This part takes a very long time
    ## Get log likelihood ratio,this is one by term maybe use in odds ratio!!!
    #set_sum_contrasts()
    #list_anova <- mclapply(list_mod,function(x){Anova(x,type="III")}, mc.preschedule=FALSE,mc.cores = num_cores) #this value is about 1/2 of deviance
    #names(list_anova) <- names_mod
    #LR_Chi_values <- unlist(lapply( list_anova ,function(x){x$`LR Chisq`})) #this value is about 1/2 of deviance
    #LR_Chi_p_values <- unlist(lapply( list_anova ,function(x){x$`Pr(>Chisq)`})) #this value is about 1/2 of deviance
    #name_list_anova <- file.path(".",paste("list_anova_",region_name[i],"_",out_suffix_s,".RData",sep=""))
    #save(list_anova,file= name_list_anova )
    
    #test<-Anova(mod,type="III")
    #test$`LR Chisq`
    
    df_val <- data.frame(mod = names_mod,
                         AIC=AIC_values,
                         deviance=res_deviance_values,
                         loglikelihood=loglikelihood_values,
                         #LR = LR_Chi_values,
                         #LR_p=LR_Chi_p_values
                         n=n_obs)
    df_val$region <-region_name[i]
    #df_val$mod <- rownames(df_val)
    rownames(df_val) <- NULL
    list_df_val[[i]] <- df_val
  }
  
  df_all <-do.call(rbind,list_df_val)
  out_filename <- file.path(out_dir,paste("tb_models_summary_accuracy_",out_suffix,".txt",sep=""))
  write.table(df_all,file=out_filename,sep=",",row.names = F)
  
  return(df_all) 
}

############### END OF SCRIPT ###################

#http://www.r-bloggers.com/r-code-example-for-neural-networks/
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
#http://data.princeton.edu/wws509/r/c6s2.html

