library(data.table)
library(tidyverse)
library(mvtnorm)
library(parallel)
library(bfast)
library(MCMCpack)
library(reticulate)
source("ccdc.R")
source("priors_myanmar.R")
library(BOCPD)
# Before starting:
## Each function will output its result to the Global Env, regardless of your
## choice for the "saveRdata" argument.
### The argument "priorRdata" asks if you want to load a saved Rdata file from
### the previous chronological script.


n.cores=16


# choose analyses to conduct
ft <- "myanmar"


# variables needed ####
set.seed(42)
country <- "myanmar" #all scripts
dataFolder <- paste0("data/", country, "/") #all scripts
inputFolder <- "geeAll" #geeWindow or geeNewDec20
#inputFolder <- "geeNewDec20"
#trainingPointsFile <- "training_points_window.csv" 
trainingPointsFile <- "training_points_all.csv"#"training_points_all.csv"
index <- "ndvi" ## eventually to include more
method = "har" ## method = "ccdc" or "har" [nonlinear harmonic]
desp = "median" ## desp = "median", "Frost", "RLF", "GMfull"
pol = "vv" ## eventually include "vh" or "ratio"
look = "ascending" ## or "descending"

#load pointIDs
base <- fread(paste0(dataFolder, trainingPointsFile))
pointids <- sort(unique(base[,pointid]))

#0. Prep S1 data ####
# MUST RUN THIS IF NEW Sentinel 1 GEE DATA
# this consolidates all the s1 speckle filters to one table
# file created = "sentinel1.csv"

# source("scripts/0_preps1.R")
# preps1_despeckle("myanmar", inputFolder=inputFolder)

#1. cloud-masking and satellite data processing ####
## variable created = onboard (list of 2 - some quick stats, and data for next steps)
## file created = "sat_data_[country].Rdata"
#source("scripts/1_cloudmasking.R")
# F1_onboard <- onboard_cloud(country, inputFolder=inputFolder, 
#                             pointData=trainingPointsFile,
#                             saveRdata=FALSE)
pointData <- trainingPointsFile


fn <- paste0(dataFolder,inputFolder,"/landsat8sr.csv")
df <- read.table(fn,sep=",",header=T)
df$date <- as.Date(df$date,format="%Y-%m-%d")

bandNames = c("SR_B4", "SR_B5","SR_B7", "QA_PIXEL", "QA_RADSAT", "SR_QA_AEROSOL")
indices <- filterThenIndex(df,sat="landsat8sr",bandNames=bandNames)

colnames(indices) <- c("valid","ndvi","swir2")

df <- data.table(df,indices)

#bring in dates to be used for model runs
base <- fread(paste0(dataFolder, pointData))
base <- base[,.(pointid, datePre, dateDist)]
setnames(base, old=c("datePre", "dateDist"),
         new=c("date_model", "date_dist"))

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

#2. calculate NDVI and SWIR
# df <- F1_onboard$sat_list$l8sr
# df <- df %>% mutate(ndvi = (B5-B4)/(B5+B4))
# df <- df %>% mutate(swir2 = B7/1e4)
# split into point ids for pure dataset
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

#df <- df %>% filter(valid_qa&valid_radsatqa&valid_aerosol&valid_cloudAndRadqa)

#4. split into pointids and sort
dfl <- split(df,df$pointid)
dfl <- lapply(dfl,function(x){
  xi <- mutate(x,t=as.numeric(x$date-min(x$date)))
  xi <- arrange(xi,date)
  mutate(xi,ti=1:nrow(xi))
  # filter
})
#5. Get priors


sample_data_list_pure <- mclapply(dfl_pure,function(dfi){
  t <- dfi$t
  period=365
  
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),(t-min(t))/(max(t)-min(t)))
  
  Y <- cbind(dfi$ndvi,dfi$swir2)
  
  find_scale(Y,X,hls_sorted=dfi$date,change_date_p = dfi$date_dist[1],ptype=3)
},mc.cores=n.cores)

priors <- get_priorsm(data_list=sample_data_list_pure,d=2,k=4)
#save(priors,file="myanmar_priors.Rdata")

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

# general priors - combine all scenario info
period=365

set.seed(53226)

sim_settings <- matrix(c(1,priors$nu,0.5,0.5,1,3),nrow=13,ncol=6,byrow=T)
sim_settings[2,1] <- 10
sim_settings[3,1] <- 0.1
sim_settings[4,1] <- 0.01
sim_settings[5,2] <- 3.01
sim_settings[6,2] <- 10
sim_settings[7,3] <- 0.9
sim_settings[8,4] <- 0.1
sim_settings[9,4] <- 0.9
sim_settings[10,5] <- 0.1
sim_settings[11,5] <- 2
sim_settings[12,6] <- 5
sim_settings[13,6] <- 1




#plan(multisession)

#cl <- makeCluster(detectCores()-1)
#registerDoParallel(cl)

#cl <- parallel::makeCluster(n_cores)

#doParallel::registerDoParallel(cl)

for(rw in 1:nrow(sim_settings)){
  lam <- sim_settings[rw,1]
  nu <- sim_settings[rw,2]
  alpha <- sim_settings[rw,3]
  pc <- sim_settings[rw,4]
  vs <- sim_settings[rw,5]
  L <- sim_settings[rw,6]
  
  ####
  piMine <- priors
  piMine$Lambda <- lam*priors$Lambda
  piMine$V <- vs*((nu-2-1)/(priors$nu-2-1))*priors$V
  piMine$nu <- nu
  
  
  
  #all_results <- foreach(n = 1:maxiters, .packages = c("mvtnorm","bfast","MCMCpack","reticulate"), .errorhandling="pass") %dopar% {
  all_results <- mclapply(dfl,function(dfi){
    # initialize results list
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
    
    
    
    
    
    # analyze using BOCPD
    start=proc.time()
    bocpd_mod <- bocpd(datapts = Yi,
                       covariates = X,
                       BOCPD = NULL,
                       par_inits = piMine,
                       Lsearch = 15,
                       Lwindow = 8,
                       Lgroup=L,
                       lambda=100,
                       cpthresh = 0.8,
                       truncRmin = 300,
                       cptimemin = 12,
                       Lm = 15,
                       cp_delay = 3,
                       kt = 1,
                       alpha=alpha,
                       pc=pc,
                       getR = FALSE,
                       getOutliers = TRUE,
                       getModels = FALSE)

    output$bocpd_time = (proc.time()-start)[3]
    
    cp_info <- extract_cps(bocpd_mod,truecp, 0.8,cpmin=12)
    cp_info
    output$bocpd_cps=cp_info$cp_dates
    output$bocpd_det=cp_info$cp_detected
    output$bocpd_lat=cp_info$cp_lat
    output$bocpd_outliers=bocpd_mod$outliers
    
    return(output)
  },mc.cores=n.cores,mc.preschedule=F)
  #}
  save(file = paste("sensitivity",rw,ft,".RData",sep=""),list=c("all_results","dfl"))
  
}
