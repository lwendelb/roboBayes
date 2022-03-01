library(mvtnorm)
library(parallel)
library(bfast)
library(MCMCpack)
library(reticulate)
#library(doParallel)
#library(foreach)
#library(future.apply)
#library(matrixNormal)
library(BOCPD)
source("ccdc.R")
source("sim_study_helpers.R")


maxiters <- 1000

n.cores=16

scens <- 1:9
#scens <- 5:8
#scens <- 9

# which analysis to collect

ft <- "sim_v2"
run.BOCPDOD <- T
run.BOCPD <- T
run.CCDC <- T
run.BFAST <- F
run.BOCPDMS <- F
run.rBOCPDMS <- F

# general priors - combine all scenario info
period=365


#t <- seq(1,20*365,length.out=20*5)
#X <- cbind(1,sin(2*pi*t/period)+cos(2*pi*t/period),t/1000)
if(all(scens==1:4)){
  set.seed(3269275)
}else{
  set.seed(45814)
}

#plan(multisession)

#cl <- makeCluster(detectCores()-1)
#registerDoParallel(cl)

#cl <- parallel::makeCluster(n_cores)

#doParallel::registerDoParallel(cl)

for(scen in scens){

####

  
#all_results <- foreach(n = 1:maxiters, .packages = c("mvtnorm","bfast","MCMCpack","reticulate"), .errorhandling="pass") %dopar% {
all_results <- mclapply(1:maxiters,function(iter,scenario=scen){
  # initialize results list
  output <- list()

  # generate data
  yrs <- 15
  npts=18
  t <- seq(1,yrs*365,length.out=yrs*npts)
  period=365
  
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),t/1000)
  
  make_scenario <-  generate_scenario(scenario,yrs=yrs,npts=npts,l1=10,l2=5)
  
  Y <- make_scenario$Y
  # priors_specific <- make_scenario$priors
  if(scenario <=4){
    priors_general <- get_priors(X,mu=0.5,seasonal1=0,seasonal2=0,rho=0.9)
  } 
  if((scenario > 4)&(scenario <=8)){
    priors_general <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
  }
  if(scenario ==9){
    priors_general <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0)
  }
  priors <- priors_general
  if(any(scenario==c(1,2,5,6))){
    newV <- priors$V
    newV[1,2] <- 0
    newV[2,1] <- 0
    priors$V <- newV
  }
  piMine <- priors
  #piMine$Lambda <- 0.01*priors$Lambda
  # 
  #priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
  
  if(run.CCDC){
    # analyze using CCDC
    start = proc.time()
    Y_ccdc <- list(Y[,1],Y[,2])
    ccdc_results <- runCCDC(Y_ccdc,t,X,ts=12)
    output$ccdc_cps = ccdc_results$cps[-1]
    output$ccdc_time = (proc.time()-start)[3]
  }
  
  if(run.BFAST){
    # analyze using BFAST
    start = proc.time()
    bfast_results <- runBFAST(Y,t,tfreq=npts,start_cp=c(6,1))
    output$bfast_cps = bfast_results$cps
    output$bfast_detected = bfast_results$cp_det
    output$bfast_time = (proc.time()-start)[3]
  }

  
  if(run.BOCPDOD){
    # analyze using BOCPD
    start=proc.time()
    
    bocpd_mod <- bocpd(datapts = Y,
                       covariates = X,
                       BOCPD = NULL,
                       par_inits = piMine,
                       Lsearch = 15,
                       Lwindow = 20,
                       Lgroup=5,
                       lambda=270,
                       cpthresh = 0.8,
                       truncRmin = 300,
                       cptimemin = 50,
                       Lm = 6,
                       cp_delay = 4,
                       kt = 1,
                       alpha=0.5,
                       pc=0.1,
                       getR = FALSE,
                       getOutliers = TRUE,
                       getModels = FALSE)

    output$bocpd_time = (proc.time()-start)[3]
    
    cp_info <- extract_cps(bocpd_mod, 0.8,cpmin=50)
    cp_info
    
    output$bocpd_cps=cp_info$cp_dates
    output$bocpd_det=cp_info$cp_detected
    output$bocpd_lat=cp_info$cp_lat
    output$bocpd_outliers=bocpd_mod$outliers
  }
  if(run.BOCPD){
    start=proc.time()
    bocpd_mod_no <- bocpd_mod <- bocpd(datapts = Y,
                                       covariates = X,
                                       BOCPD = NULL,
                                       par_inits = piMine,
                                       Lsearch = 15,
                                       Lwindow = 20,
                                       Lgroup=5,
                                       lambda=270,
                                       cpthresh = 0.8,
                                       truncRmin = 300,
                                       cptimemin = 50,
                                       Lm = 7,
                                       cp_delay = 4,
                                       kt = 1,
                                       alpha=0.5,
                                       pc=0.5,
                                       getR = FALSE,
                                       getOutliers = FALSE,
                                       getModels = FALSE)
    
    
    output$bocpd2_time = (proc.time()-start)[3]
    
    cp_info_no <- extract_cps(bocpd_mod_no, 0.8,cpmin=50)
    
    output$bocpd2_cps=cp_info_no$cp_dates
    output$bocpd2_det=cp_info_no$cp_detected
    output$bocpd2_lat=cp_info_no$cp_lat
  }
  
  # analyze with robust bocpd
  #setwd("C:/Users/Laura/Documents/GitHub/rbocpdms/rbocpdms-master")
  
  #write.table(Y,"C:/Users/Laura/Documents/GitHub/rbocpdms/rbocpdms-master/Data/well log/sim.txt",row.names=F,col.names=F,sep=",")
  #setwd("C:/Users/Laura/SkyDrive/Documents/NCSU/Spring 2021/CPD/R")
  
  #write.table(Y,"sim.txt",row.names=F,col.names=F,sep=",")
  
  if(run.rBOCPDMS){
    assign("Y", Y, envir = globalenv())
    start=proc.time()
    py_run_file("sim_study_rbocpdms.py")
    output$rbocpdms_time = (proc.time()-start)[3]
    
    rbocpdms_cps <- lapply(1:269,function(i){sapply(py$detector$CPs[[i]],function(x){
      x[[1]][1]
    })})
    
    output$rbocpdms <- rbocpdms_cps
  }
  if(run.BOCPDMS){
    assign("Y", Y, envir = globalenv())
    start=proc.time()
    py_run_file("sim_study_rbocpdms.py")
    output$bocpdms_time = (proc.time()-start)[3]
    
    bocpdms_cps <- lapply(1:269,function(i){sapply(py$detector$CPs[[i]],function(x){
      x[[1]][1]
    })})
    bocpdms_cps[[269]]
    output$bocpdms <- bocpdms_cps
  }
  
  return(output)
},mc.cores=n.cores,mc.preschedule=F)
#}
save(file = paste("ss",scen,ft,".RData",sep=""),list=c("all_results"))

}

#stopCluster(cl)

# par(mfrow=c(2,1),mar=c(2,4,4,2))
# #plot(Y[,1])
# #plot(Y[,2])
# 
# plot(Y[,1],ylim=c(-0.2,1),xlab="Time point",ylab="Signal 1",main="Detection of simulated correlation change at 181 using both signals")
# abline(v=unique(unlist(bocpd_mod$cp_inds))[-1],col="blue")
# abline(v=bocpd_mod$outliers,col="red")
# legend("bottomleft",legend=c("Multivariate analysis","Univariate Signal 1 analysis"),
#        lty=c(1,2),col=c("blue","red"))
# plot(Y[,2],ylim=c(-0.2,1),xlab="Time point",ylab="Signal 2")
# abline(v=unique(unlist(bocpd_mod$cp_inds))[-1],col="blue")
# abline(v=bocpd_mod$outliers,col="red")
# legend("bottomleft",legend=c("Multivariate analysis","Univariate Signal 2 analysis"),
#        lty=c(1,2),col=c("blue","red"))



# make_scenario <-  generate_scenario(scenario,yrs=yrs,npts=npts,l1=10,l2=5)
# 
# Y <- make_scenario$Y
# priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
# 
# bocpd_mod <- runBOCPD(datapts=Y,
#                       covariates=X,
#                       hazard=function(r){geom_hazard(r,1000)},
#                       par_inits=list(B=priors$B,
#                                      V=priors$V,
#                                      nu=priors$nu,
#                                      Lambda=priors$Lambda,
#                                      p=0, ps=0, d = 2,k=4,kt=1),
#                       truncRthresh=1e-4,truncRmin=300,
#                       cpthresh=0.8,cptimemin=50,Lgroup=5,Lsearch=15,
#                       Lwindow=20,cp_delay=3,verbose=F,
#                       getR=T,getMean=F,getR_m=T)
# cp_info <- extract_cps(bocpd_mod, 0.8,cpmin=50)
# cp_info
# bocpd_mod$outliers
# abline(v=184)
