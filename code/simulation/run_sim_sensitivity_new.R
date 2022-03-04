library(mvtnorm)
library(parallel)
library(bfast)
library(MCMCpack)
library(reticulate)
library(roboBayes)
source("ccdc.R")
source("sim_study_helpers.R")


maxiters <- 1000

n.cores=16

scen<- 7

# which analysis to collect

ft <- "sensitivity"

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


for(rw in 1:nrow(sim_settings)){
  mu <- sim_settings[rw,1]
  lam <- sim_settings[rw,2]
  nu <- sim_settings[rw,3]
  alpha <- sim_settings[rw,4]
  pc <- sim_settings[rw,5]
  L <- sim_settings[rw,6]
  
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
      # analyze using roboBayes
      start=proc.time()
      
      bocpd_mod <- roboBayes(datapts = Y,
                         covariates = X,
                         roboBayes = NULL,
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
  save(file = paste("results/ss",rw,ft,".RData",sep=""),list=c("all_results"))
  
}

