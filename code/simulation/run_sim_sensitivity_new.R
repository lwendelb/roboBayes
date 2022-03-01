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

#scens <- 1:4
#scens <- 5:8
scen<- 7

# which analysis to collect

ft <- "sensitivity_v2"

# general priors - combine all scenario info
period=365

set.seed(53226)

sim_settings <- matrix(c(0.5,1,20,0.5,0.1,5),nrow=11,ncol=6,byrow=T)
sim_settings[2,1] <- 0.4
sim_settings[3,2] <- 0.01
sim_settings[4,2] <- 100
sim_settings[5,3] <- 4
sim_settings[6,3] <- 50
sim_settings[7,4] <- 0.9
sim_settings[8,5] <- 0.5
sim_settings[9,5] <- 0.9
sim_settings[10,6] <- 1
sim_settings[11,6] <- 3




#plan(multisession)

#cl <- makeCluster(detectCores()-1)
#registerDoParallel(cl)

#cl <- parallel::makeCluster(n_cores)

#doParallel::registerDoParallel(cl)

for(rw in 1:nrow(sim_settings)){
  mu <- sim_settings[rw,1]
  lam <- sim_settings[rw,2]
  nu <- sim_settings[rw,3]
  alpha <- sim_settings[rw,4]
  pc <- sim_settings[rw,5]
  L <- sim_settings[rw,6]
  
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
      priors_general <- get_priors(X,mu=mu,seasonal1=0,seasonal2=0,rho=0.9)
    } 
    if((scenario > 4)&(scenario <=8)){
      priors_general <- get_priors(X,mu=mu,seasonal1=0.1,seasonal2=0.04,rho=0.9)
    }
    if(scenario ==9){
      priors_general <- get_priors(X,mu=mu,seasonal1=0.1,seasonal2=0.04,rho=0)
    }
    priors <- priors_general

    ####
    piMine <- priors
    piMine$Lambda <- lam*priors$Lambda
    piMine$V <- ((nu-2-1)/(priors$nu-2-1))*priors$V
    piMine$nu <- nu
      # analyze using BOCPD
      start=proc.time()
      
      bocpd_mod <- bocpd(datapts = Y,
                         covariates = X,
                         BOCPD = NULL,
                         par_inits = piMine,
                         Lsearch = 15,
                         Lwindow = 20,
                         Lgroup=L,
                         lambda=270,
                         cpthresh = 0.8,
                         truncRmin = 300,
                         cptimemin = 50,
                         Lm = 6,
                         cp_delay = 4,
                         kt = 1,
                         alpha=alpha,
                         pc=pc,
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
    
    return(output)
  },mc.cores=n.cores,mc.preschedule=F)
  #}
  save(file = paste("ss",rw,ft,".RData",sep=""),list=c("all_results"))
  
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
