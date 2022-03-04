library(data.table)
library(tidyverse)
library(mvtnorm)
library(parallel)
library(MCMCpack)
library(reticulate)
source("ccdc.R")
source("priors_myanmar.R")
library(roboBayes)
# Before starting:


n.cores=16


# choose analyses to conduct
ft <- "myanmar"
run.roboBayes <- T
run.BOCPD <- T
run.CCDC <- T
run.BOCPDMS <- F
run.rBOCPDMS <- F


# variables needed ####
set.seed(42)
country <- "myanmar" #all scripts
dataFolder <- paste0("data/", country, "/") #all scripts
inputFolder <- "analysis_points"
trainingPointsFile <- "analysis_points.csv"

#load pointIDs
base <- fread(paste0(dataFolder, trainingPointsFile))
pointids <- sort(unique(base[,pointid]))

pointData <- trainingPointsFile

###### Data prep
fn <- paste0(dataFolder,inputFolder,"/landsat8sr.csv")
df <- read.table(fn,sep=",",header=T)
df$date <- as.Date(df$date,format="%Y-%m-%d")

bandNames = c("SR_B4", "SR_B5","SR_B7", "QA_PIXEL", "QA_RADSAT", "SR_QA_AEROSOL")
indices <- filterThenIndex(df,sat="landsat8sr",bandNames=bandNames)

colnames(indices) <- c("valid","ndvi","swir2")

df <- data.table(df,indices)

#bring in dates to be used for model runs
base <- fread(paste0(dataFolder, pointData))
base <- base[,.(pointid, datePre, dateDist,distPercL8)]
setnames(base, old=c("datePre", "dateDist","distPercL8"),
         new=c("date_model", "date_dist","dist_perc"))

ids <- sort(unique(base[,pointid]))

#update dates
date_fn <- function(dae){
  date_true <- ifelse(is.na(as.Date(dae)) | (Sys.Date() - as.Date(dae)) > 10000, 
                      as.Date(dae, format="%d/%m/%Y"),
                      as.Date(dae))
  date_true <- as.Date(date_true, origin="1970-01-01")
  return(date_true)
}

base <- base[, `:=` (date_model = date_fn(date_model),
                     date_dist = date_fn(date_dist))]

get_moddate <- function(sat_list=sat_list){
  dat <- sat_list
  setkey(dat, pointid)
  setkey(base, pointid)
  dat <- dat[base]
}

df <- get_moddate(data.table(df))

# split into point ids for pure dataset (no corrupted points)
df_pure <- df %>% mutate(goodquality=valid)
#4. filter out bad qa
df_pure <- df_pure %>% filter(goodquality)

dfl_pure <- split(df_pure,df_pure$pointid)
dfl_pure <- lapply(dfl_pure,function(x){
  xi <- mutate(x,t=as.numeric(x$date-min(x$date)))
  xi <- arrange(xi,date)
  mutate(xi,ti=1:nrow(xi))
  # filter
})
#3. Keep in some bad points
corrupted <- sample(which(!df$valid&!is.na(df$ndvi)&!is.na(df$swir2)),200)
df <- df %>% mutate(goodquality=valid)
df$goodquality[corrupted] <- T
#4. filter out bad qa
df <- df %>% filter(goodquality)

#4. split into pointids and sort
dfl <- split(df,df$pointid)
dfl <- lapply(dfl,function(x){
  xi <- mutate(x,t=as.numeric(x$date-min(x$date)))
  xi <- arrange(xi,date)
  mutate(xi,ti=1:nrow(xi))
  # filter
})
#5. Get priors

# fit models on pure data to estimate prior hyperparameters
sample_data_list_pure <- mclapply(dfl_pure,function(dfi){
  t <- dfi$t
  period=365
  
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),(t-min(t))/(max(t)-min(t)))
  
  Y <- cbind(dfi$ndvi,dfi$swir2)
  
  find_scale(Y,X,hls_sorted=dfi$date,change_date_p = dfi$date_dist[1],ptype=3)
},mc.cores=n.cores)

priors <- get_priorsm(data_list=sample_data_list_pure,d=2,k=4)

piMine <- priors

# fit models and residuals on corrupted data to get residuals
sample_data_list <- mclapply(dfl,function(dfi){
  t <- dfi$t
  period=365
  
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),(t-min(t))/(max(t)-min(t)))
  
  Y <- cbind(dfi$ndvi,dfi$swir2)

  
  find_scale(Y,X,hls_sorted=dfi$date,change_date_p = dfi$date_dist[1],ptype=3)
},mc.cores=n.cores)

# add residuals into dfls
dfl <- mapply(function(dfi,x){
  dfi$e1 <- x$prior_vals$ehat[,1]
  dfi$e2 <- x$prior_vals$ehat[,2]
  return(dfi)
},dfi=dfl,x=sample_data_list,SIMPLIFY = F)

all_results <- mclapply(dfl,function(dfi){
  output <- list()
  # make covariate matrix
  t <- dfi$t
  period=365
  n <- length(t)
  
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),(t-min(t))/(max(t)-min(t)))
  
  Yi <- cbind(dfi$ndvi,dfi$swir2)
  Y <- cbind(dfi$e1,dfi$e2)
  
  # pull out change point number
  truecp <- dfi$ti[which(dfi$date_dist<=dfi$date)[1]]
  
  if(run.CCDC){
    # analyze using CCDC
    start = proc.time()
    Y_ccdc <- list(Yi[,1],Yi[,2])
    ccdc_results <- runCCDC(Y_ccdc,t,X,ts=12)
    output$ccdc_cps = ccdc_results$cps[-1]
    output$ccdc_time = (proc.time()-start)[3]
  }
  
  if(run.roboBayes){
    # analyze using BOCPD
    start=proc.time()
    bocpd_mod <- roboBayes(datapts = Yi,
                          covariates = X,
                          RoboBayes = NULL,
                          par_inits = piMine,
                          Lsearch = 15,
                          Lwindow = 8,
                          Lgroup=3,
                          lambda=100,
                          cpthresh = 0.8,
                          truncRmin = 300,
                          cptimemin = 12,
                          Lm = 8,
                          cp_delay = 3,
                          kt = 1,
                          pc=0.5,
                          alpha=0.5,
                          getR = FALSE,
                          getOutliers = TRUE,
                          getModels = FALSE)
    output$bocpd_time = (proc.time()-start)[3]
    
    cp_info <- extract_cps(bocpd_mod,truecp, 0.8,cpmin=12)
    output$bocpd_cps=cp_info$cp_dates
    output$bocpd_det=cp_info$cp_detected
    output$bocpd_lat=cp_info$cp_lat
    output$bocpd_outliers=bocpd_mod$outliers
  }
  if(run.BOCPD){
    start=proc.time()
    bocpd_mod_no <- roboBayes(datapts = Yi,
                   covariates = X,
                   RoboBayes = NULL,
                   par_inits = piMine,
                   Lsearch = 15,
                   Lwindow = 8,
                   Lgroup=3,
                   lambda=100,
                   cpthresh = 0.8,
                   truncRmin = 300,
                   cptimemin = 12,
                   Lm = 4,
                   cp_delay = 3,
                   kt = 1,
                   getR = FALSE,
                   getOutliers = FALSE,
                   getModels = FALSE)
    output$bocpd2_time = (proc.time()-start)[3]
    
    
    cp_info_no <- extract_cps(bocpd_mod_no,truecp, 0.8,cpmin=12)
    
    output$bocpd2_cps=cp_info_no$cp_dates
    output$bocpd2_det=cp_info_no$cp_detected
    output$bocpd2_lat=cp_info_no$cp_lat
  }
  
  # analyze with robust bocpd
  if(run.rBOCPDMS){
    assign("Y", Y, envir = globalenv())
    start=proc.time()
    py_run_file("myanmar_rbocpdms.py")
    output$rbocpdms_time = (proc.time()-start)[3]
    
    rbocpdms_cps <- lapply(1:(n-1),function(i){sapply(py$detector$CPs[[i]],function(x){
      x[[1]][1]
    })})
    
    output$rbocpdms <- rbocpdms_cps
  }
  if(run.BOCPDMS){
    assign("Y", Y, envir = globalenv())
    start=proc.time()
    py_run_file("myanmar_bocpdms.py")
    output$bocpdms_time = (proc.time()-start)[3]
    
    bocpdms_cps <- lapply(1:(n-1),function(i){sapply(py$detector$CPs[[i]],function(x){
      x[[1]][1]
    })})
    output$bocpdms <- bocpdms_cps
  }
  
  return(output)
},mc.cores=n.cores,mc.preschedule=F)
#}
save(file = paste("results/",ft,".RData",sep=""),list=c("all_results","dfl"))

