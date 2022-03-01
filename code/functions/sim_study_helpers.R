# data generation function
generate <- function(mu,rho,t,
                     seasonal=F,noerror=F,outlier=F,period=365){

  
  d=2
  k=4
  # B0 <- matrix(c(1,0.1,0.2,0,0.5,-0.1,-0.2,0),4,2)
  Lambda0 <- 0.001*diag(c(100,10000,10000,10000))
  nu0 <- 20
  V0 <- (nu0-d-1)*matrix(0.001*c(1,rho,rho,1),2,2)
  
  X <- cbind(1,
             sin(2*pi*t/period),
             cos(2*pi*t/period),
             t/1000)
  Sigma <- riwish(nu0,V0)
  #Beta <- rmatnorm(s=1,M=B0,U=solve(Lambda0),V=Sigma)
  if(seasonal){
    B0 <- rbind(c(mu),matrix(c(.1,.04,0,.1,.04,0),3,2))
    Beta <- matrix(rmvnorm(1,c(B0),kronecker(Sigma,solve(Lambda0))),k,d)
    
    Beta[1,] <- mu
  }
  else{
    B0 <- rbind(c(mu),matrix(0,3,2))
    Beta <- B0
  }

  
  E <- rmvnorm(length(t),rep(0,2),Sigma)
  Y <- X%*%Beta + E
  
  #plot(t,Y[,1])
  if(noerror){
    return(Y-E)
  }else{
    return(Y)
  }
  
}

get_priors <- function(X,mu=0.5,seasonal1=0.025,seasonal2=-0.01,rho=0,nu=10){
  B0 <- matrix(c(rep(mu,2),rep(seasonal1,2),rep(seasonal2,2),0,0),4,2,byrow=T)
  
  d=2
  k=4
  # B0 <- matrix(c(1,0.1,0.2,0,0.5,-0.1,-0.2,0),4,2)
  Lambda0 <- 0.001*diag(c(100,10000,10000,10000))
  nu0 <- 20
  V0 <- (nu0-d-1)*matrix(0.001*c(1,rho,rho,1),2,2)
  
  priors <- list(B=B0,
                 Lambda=Lambda0,
                 V=V0,
                 nu = nu0)
}

priors_sensitivity <- function(i){
  if(i==1){
    return(list())
  }
}

generate_scenario <- function(scenario,yrs=10,npts=20,noerror=F,nu=10,l1=10,l2=NULL){
  if(is.null(l2)){
    l2 <- yrs-l1
  }
  if(is.null(yrs)){
    yrs <- l1+l2
  }
  t <- seq(1,yrs*365,length.out=npts*yrs)
  X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),t/1000)
  t1 <- t[1:(npts*l1)]
  t2 <- t[(npts*l1+1):(npts*(l1+l2))]
  
  outlier_val <- c(0.8,0.1)
  if(scenario==1){
    y1 <- generate(0.5,0,t1,seasonal=F,noerror=noerror)
    y2 <- generate(0.4,0,t2,seasonal=F,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0,seasonal2=0,rho=0)
  }
  if(scenario==2){
    y1 <- generate(0.5,0,t1,seasonal=F,noerror=noerror)
    y2 <- generate(0.3,0,t2,seasonal=F,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0,seasonal2=0,rho=0)
  }
  if(scenario==3){
    y1 <- generate(0.5,0.9,t1,seasonal=F,noerror=noerror)
    y2 <- generate(0.4,0.9,t2,seasonal=F,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0,seasonal2=0,rho=0.9)

  }
  if(scenario==4){
    y1 <- generate(0.5,0.9,t1,seasonal=F,noerror=noerror)
    y2 <- generate(0.3,0.9,t2,seasonal=F,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0,seasonal2=0,rho=0.9)

  }
  if(scenario==5){
    y1 <- generate(0.5,0,t1,seasonal=T,noerror=noerror)
    y2 <- generate(0.4,0,t2,seasonal=T,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0)
    
  }
  if(scenario==6){
    y1 <- generate(0.5,0,t1,seasonal=T,noerror=noerror)
    y2 <- generate(0.3,0,t2,seasonal=T,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0)
    
  }
  if(scenario==7){
    y1 <- generate(0.5,0.9,t1,seasonal=T,noerror=noerror)
    y2 <- generate(0.4,0.9,t2,seasonal=T,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
    
  }
  if(scenario==8){
    y1 <- generate(0.5,0.9,t1,seasonal=T,noerror=noerror)
    y2 <- generate(0.3,0.9,t2,seasonal=T,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
    
  }
  if(scenario==9){
    y1 <- generate(0.5,-0.5,t1,seasonal=T,noerror=noerror)
    y2 <- generate(0.5,0.5,t2,seasonal=T,noerror=noerror)
    Y <- rbind(y1,y2)
    oi <- sample(round(l1*npts/2+1):length(t),1)
    Y[oi,] <- outlier_val
    priors <- get_priors(X,mu=0.5,seasonal1=0.1,seasonal2=0.04,rho=0.9)
  }
  return(list(Y=Y,priors=priors,outlier=oi))
}

extract_cps <- function(bocpd_mod,cpthresh,cpmin){
  
  
  cps_unique <- lapply(1:length(cpthresh),function(x){
    unique(unlist(bocpd_mod$cp_inds[[x]]))[-1]
  })
  cps_detected <- lapply(1:length(cpthresh),function(j){
    cpd <- rep(0,length(cps_unique[[j]]))
    if(length(cps_unique[[j]])>0){
      for(i in 1:length(cps_unique[[j]])){
        cpd[i] <- cpmin + which(unlist(lapply(bocpd_mod$cp_inds[[j]][-1],function(x){cps_unique[[j]][i]==x})))[1]
      }
    }
    return(cpd)
  })

  
  cp_lat <- lapply(1:length(cpthresh),function(x){
    cps_detected[[x]]-cps_unique[[x]]})
  
  
  # 1st column is the change point, 2nd column is the initial 
  # date of detection for the changepoint
  
  return(list(cp_dates=cps_unique,cp_detected=cps_detected,cp_lat=cp_lat))
}


extract_cps <- function(bocpd_mod,truecp,cpthresh,cpmin){
  # first get cps
  cps <- unique(c(bocpd_mod$cpInds))[-1]
  # eliminate ones in the training period
  if(any(cps<cpmin)){
    cps <- cps[cps>=cpmin]
  }
  # get latency
  cpd <- 0*cps
  if(length(cps)>0){
    for(cpi in 1:length(cps)){
      cpd[cpi] <- bocpd_mod$time-1-length(bocpd_mod$cpInds)+which(sapply(bocpd_mod$cpInds,function(xi){
        any(cps[cpi]==xi)
      }))[1]
    }
  }
  
  cp_lat <- cpd-cps
  # 1st column is the change point, 2nd column is the initial 
  # date of detection for the changepoint
  
  return(list(cp_dates=cps,cp_detected=cpd,cp_lat=cp_lat))
}

# extract_cps <- function(bocpd_mod,truecp,cpthresh,cpmin){
#   # first get cps
#   cps <- unique(unlist(bocpd_mod$cp_inds))[-1]
#   # eliminate ones in the training period
#   if(any(cps<cpmin)){
#     cps <- cps[cps>=cpmin]
#   }
#   # get latency
#   cpd <- 0*cps
#   if(length(cps)>0){
#     for(cpi in 1:length(cps)){
#       cpd[cpi] <- bocpd_mod$time-1-length(bocpd_mod$cp_inds[[1]])+which(sapply(bocpd_mod$cp_inds[[1]],function(xi){
#         any(cps[cpi]==xi)
#       }))[1]
#     }
#   }
#   
#   cp_lat <- cpd-cps
#   # 1st column is the change point, 2nd column is the initial 
#   # date of detection for the changepoint
#   
#   return(list(cp_dates=cps,cp_detected=cpd,cp_lat=cp_lat))
# }


