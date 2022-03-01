library(nlme)
#library(matrixNormal)
library(MCMCpack)

find_scale <- function(Y,X,hls_sorted,change_date_p,ptype){

  d = ncol(Y)
  k = 4
  
  all_t <- as.integer(hls_sorted)
  
  # eliminate unobserved timepoints
  Y_ind_obs <- 1:nrow(Y)
  Y_obs <- Y
  hls_obs <- hls_sorted[Y_ind_obs] # observed dates
  t <- as.integer(hls_obs)
  
  # extract stable period data
  
  y_stable <- as.matrix(Y_obs[hls_obs <= change_date_p,])
  date_stable <- hls_obs[hls_obs <= change_date_p]
  t_stable <- as.integer(date_stable)
  
  #ysc <- Y_obs
  
  # fit nonlinear models to gather prior info (if enough points are available)
  # fit model using linear model
  #t_stable <- (t_stable-min(all_t))
  if(nrow(y_stable)>5){
    ymn <- colMeans(y_stable)
    ysd <- apply(y_stable,2,sd)
    if(ptype==1){
      ysc <- sapply(1:ncol(Y_obs),function(i){
        (Y_obs[,i]-ymn[i])/ysd[i]
      })
    }
    if((ptype==2) | (ptype==3)){
      ysc <- sapply(1:ncol(Y_obs),function(i){
        (Y_obs[,i])#/ysd[i]
      })
    }
    x1 = 1
    x2 = sin((2*pi/365)*t_stable)
    x3 = cos((2*pi/365)*t_stable)
    x4 = (t_stable-min(all_t))/(max(all_t)-min(all_t))
    
    X_stable <- (cbind(x1,x2,x3,x4))
    
    #Vi <- solve(kronecker(diag(nrow(y_stable)),cov(y_stable)))
    
    Y_stable <- c(t(ysc[hls_obs <= change_date_p,]))
    
    X_stable2 <- do.call(rbind,lapply(1:nrow(X_stable),
                                      function(x){
                                        t(cbind(c(X_stable[x,],rep(0,k)),
                                                c(rep(0,k),X_stable[x,])))}))
    t_stable2 <- c(sapply(1:length(date_stable),function(x){rep(t_stable[x],d)}))
    
    id <- rep(seq(1,d),nrow(y_stable))
    
    dat <- data.frame(t_stable2,X_stable2,Y_stable,id)
    
    colnames(dat) <- c("ts2","x11","x12","x13","x14","x21","x22","x23","x24","ys2","id")
    
    ##### gls
    doit <- function(){
      fm1 <- (gls(ys2 ~ -1 + x11 + x12+x13+x14 + x21 +x22+x23+x24, data=dat,
                  correlation = corSymm(form = ~ 1 | ts2),
                  weights=varIdent(form = ~ 1| id),
                  method='ML',
                  control = glsControl(returnObject=TRUE,opt='nlminb')))
      # from this, we can get Bhat, Vhat
      Bhat <- matrix(fm1$coefficients,k,d)
      
      corhat <- corMatrix(fm1$modelStruct$corStruct)[[1]]
      vvec <- sqrt(diag(coef(fm1$modelStruct$varStruct, uncons = FALSE, allCoef = TRUE)))
      Sigmahat <- fm1$sigma^2*t(vvec)%*%corhat%*%vvec
      
      # get predictions at all data
      Y_vec <- c(t(ysc))
      
      x1 = 1
      x2 = sin((2*pi/365)*all_t)
      x3 = cos((2*pi/365)*all_t)
      x4 = (all_t-min(all_t))/(max(all_t)-min(all_t))
      
      X_pred <- (cbind(x1,x2,x3,x4))
      
      X_pred2 <- do.call(rbind,lapply(1:nrow(X_pred),
                                        function(x){
                                          t(cbind(c(X_pred[x,],rep(0,k)),
                                                  c(rep(0,k),X_pred[x,])))}))
      t_pred2 <- c(sapply(1:length(hls_sorted),function(x){rep(all_t[x],d)}))
      
      id <- rep(seq(1,d),nrow(Y))
      
      dat <- data.frame(t_pred2,X_pred2,Y_vec,id)
      
      colnames(dat) <- c("ts2","x11","x12","x13","x14","x21","x22","x23","x24","ys2","id")
      
      Y_pred <- predict(fm1,dat)
      
      ehat=matrix(Y_pred-Y_vec,nrow(Y),d,byrow=T)
      
      par_vals <- list(Sigmahat=Sigmahat,Bhat=Bhat,ymn=ymn,ysd=fm1$sigma*diag(vvec),ehat=ehat)
      return(par_vals)
    }
    par_vals_try <- try(
      doit(),silent=T
    )
    if(inherits(par_vals_try,"try-error")){
      ehat <- Y-matrix(colMeans(Y),nrow(Y),ncol(Y),byrow=T)
      par_vals <- list(Sigmahat=matrix(NA,d,d),Bhat=matrix(NA,k,d),
                       ymn=rep(NA,d),ysd=rep(NA,d),ehat=ehat)
    }
    else{
      par_vals <- par_vals_try
    }
    
  }else{
    ehat <- Y-matrix(colMeans(Y),nrow(Y),ncol(Y),byrow=T)
    
    par_vals <- list(Sigmahat=matrix(NA,d,d),Bhat=matrix(NA,k,d),
                     ymn=rep(NA,d),ysd=rep(NA,d),ehat=ehat)
  }
  
  return(list(prior_vals=par_vals))
}


get_priorsm <- function(data_list,k=4,d=2){
  Sigmas <- lapply(data_list,function(x){
    x$prior_vals$Sigmahat})
  Betas <- lapply(data_list,function(x){
    x$prior_vals$Bhat})
  n <- sum(sapply(Betas,function(x){
    !any(is.na(x))
  }))
  
  logdetV <- sum(sapply(Sigmas,function(x){
    log(det(x))
  }),na.rm=T)
  meanV <- matrix(rowMeans(sapply(Sigmas,function(x){
    c(x)
  }),na.rm=T),d,d)
  f_nu <- function(nu,d,n,logdetV,meanV){
    val1 <- (n/2)*log(det((nu-d-1)*meanV))
    val2 <- -n*d*log(2)/2
    val3 <- -0.5*logdetV
    val4 <- -(n/2)*sum(sapply(1:d,function(j){
      digamma(nu/2+(1-j)/2)
    }))
    val1+val2+val3+val4
  }
  
  
  doit <- function(){
    root_result <- uniroot(f_nu,interval=c(d+1.01,200),d=d,n=n,logdetV=logdetV,meanV=meanV)
    
    nuhat <- root_result$root
    
    Vhat <- meanV*(nuhat-d-1)
    B <- t(sapply(Betas,function(x){c(x)}))
    Bhat <- matrix(colMeans(B,na.rm=T),k,d)
    
    
    LambdaInvhat <- (1/d)*matrix(colMeans(t(mapply(function(Bi,Sigmai,Bhat){
      val1 <- Bi-Bhat
      val1%*%solve(Sigmai)%*%t(val1)
    },Bi=Betas,Sigmai=Sigmas,MoreArgs = list(Bhat=Bhat))),na.rm=T),k,k)
    
    Lambdahat <- solve(LambdaInvhat)
    
    return(list(B=Bhat,V=Vhat,Lambda=Lambdahat,nu=nuhat))
  }
  par_vals_try <- try(
    doit(),silent=T
  )
  if(inherits(par_vals_try,"try-error")){
    par_vals <- list(B=matrix(NA,k,d),V=matrix(NA,d,d),
                     Lambda=matrix(NA,k,k),nu=NA)
  } else{
    par_vals <- par_vals_try
  }
  
  return(par_vals)
}

# Functions for QA/QC filtering and calculating ndvi of matrix outputs
parseL8SR_pixel <- function(x){
  # QA_PIXEL: We only want to keep the best quality images, so only those that are
  ## clear and have low confidences
  
  # Binary
  ## Bit 0 - if pixel is fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  ## Bit 1 - if dilated cloud, then true
  dilatedCloud <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  ## Bit 2 - if cirrus, then true
  cirrus <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  ## Bit 3 - if cloud, then true
  cloud <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  ## Bit 4 - if cloud shadow, then true
  cloudShadow <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  ## Bit 5 - if snow, then true
  snow <- ifelse(bitwAnd(bitwShiftR(x,5),1), TRUE, FALSE)
  ## Bit 6 - if clear, then true
  clear <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  ## Bit 7 - if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x,7),1), TRUE, FALSE)
  
  # Confidences
  ## Confidences should be interpreted as the answer to the question: "What are the chances
  ### I will see X outside?", with X being cloud, cloud shadow, etc. 
  
  ## Bits 8-9 - if cloud conf low or no level set, then false
  ### 0=no level set, 1=low, 2=medium, 3=high
  cloudConf <- ifelse(bitwAnd(bitwShiftR(x,8), 3) == 1, FALSE, TRUE)
  ## Bits 10-11 - if cloud shadow confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cloudShadowConf <- ifelse(bitwAnd(bitwShiftR(x,10), 3) == 1, FALSE, TRUE)
  ## Bits 12-13 - if snow/ice confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  snowConf <- ifelse(bitwAnd(bitwShiftR(x, 12), 3) == 1, FALSE, TRUE)
  ## Bits 14-15 - if low cirrus confidence, then false; 
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cirrusConf <- ifelse(bitwAnd(bitwShiftR(x,14), 3) == 1, FALSE, TRUE)
  
  return(list(fill=fill, dilatedCloud=dilatedCloud, cirrus=cirrus, cloud=cloud, 
              cloudShadow=cloudShadow, snow=snow, clear=clear, water=water,
              cloudConf=cloudConf, cloudShadowConf=cloudShadowConf, 
              snowConf=snowConf, cirrusConf=cirrusConf)
  )
}
parseL8SR_radsat <- function(x){
  # QA_RADSAT: We only want best images, so no saturation and no occlusion
  
  #is saturated?
  b1 <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  b2 <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  b3 <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  b4 <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  b5 <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  b6 <- ifelse(bitwAnd(bitwShiftR(x,5), 1), TRUE, FALSE)
  b7 <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  #band 8 is not used
  b9 <- ifelse(bitwAnd(bitwShiftR(x,8), 1), TRUE, FALSE)
  terrainOcclusion <- ifelse(bitwAnd(bitwShiftR(x,11), 1), TRUE, FALSE)
  
  return(list(
    b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, b6=b6, b7=b7, b9=b9,
    terrainOcclusion=terrainOcclusion
  ))
}
parseL8SR_aerosol <- function(x){
  # SR_QA_AEROSOL: We want best images, so no fill, no water, and low aerosol
  ## difference with climatology if correction applied (see user guide link in
  ## L8 section below)
  
  # Bit 0: if fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  # Bit 2: if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x, 2), 1), TRUE, FALSE)
  #aerosol level; 0=climatology (no correction), 1=low, 2=med, 3=high)
  aerosolLow <- ifelse(bitwAnd(bitwShiftR(x, 6), 3) < 2, TRUE, FALSE)
  return(list(fill=fill, water=water, aerosolLow=aerosolLow))
}

qa_remap <- function(G, qa=""){
  if(qa=="scl"){
    # https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/level-2a/algorithm
    return(ifelse(G %in% c(4,5), "good", "bad"))
  }
  
  if(qa=="pixel"){
    bit_output <- parseL8SR_pixel(G)
    cond_false <- c("fill", "dilatedCloud", "cirrus", "cloud", "cloudShadow", 
                    "snow", "water", "cloudConf", "cloudShadowConf", 
                    "snowConf", "cirrusConf")
    cond_true <- c("clear")
    
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
  
  if(qa=="radsat"){
    bit_output <- parseL8SR_radsat(G)
    cond_false <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7",
                    "b9", "terrainOcclusion")
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0, "good", "bad"))
  }
  
  if(qa=="sr_aerosol"){
    ### pixel = 0 && radsat = 0. We are not filtering out based on aerosol due to low numbers of obs
    bit_output <- parseL8SR_aerosol(G)
    cond_false <- c("fill", "water")
    cond_true <- c("aerosolLow")
    #in other words, the data is valid (TRUE), if aerosol is low.
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
}

# assuming the bands are in order (e.g., Band 4 then Band 5)
calcIndex <- function(dataMat, s=sat, b=bandNames){
  
  # no scaling for S2 bc already done in L2A conversion
  bandA <- dataMat[[b[1]]]
  bandB <- dataMat[[b[2]]]
  bandC <- dataMat[[b[3]]]
  
  # scaling comes from Table 6-1 https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1619_Landsat-8-Collection2_Level-2_Science-Product-Guide-v3.pdf
  if(grepl("landsat", s)){
    bandA <- bandA*0.0000275 + (-0.2)
    bandB <- bandB*0.0000275 + (-0.2)
    bandC <- bandC*0.0000275 + (-0.2)
  }
  
  return(cbind((bandB - bandA) / (bandB + bandA),bandC))
}


# apply QA/QC to the bands individually, then calculate ndvi. This returns a matrix of ndvi
filterThenIndex <- function(input, sat=names(vrts)[X], bandNames=bands){
  qaBands <- bandNames[!grepl("B", bandNames)]
  indexBands <- bandNames[grepl("B", bandNames)]
  
  valid <- rep(T,nrow(input))
  
  if(grepl("landsat", sat)){
    #qaBands <- qaBands[!grepl("AEROSOL", qaBands)]
    qaNames <- gsub("QA_", "", qaBands)
    
    #filter out values that are not valid (outside of valid range based on user guide - see above in calcIndex)
    for(j in 1:length(indexBands)){
      valid[input[[bandNames[j]]] < 7273 | input[[bandNames[j]]] > 43636] <- F
    }
  } else {
    qaNames <- qaBands
  }
  
  for(i in 1:length(qaBands)){
    look <- data.table(v=sort(unique(as.numeric(input[[qaBands[i]]]))))
    look[, qa := sapply(v, qa_remap, qa=tolower(qaNames[i]))] # apply qa_map function
    
    valid[input[[qaBands[i]]] %in% look[qa == "bad", v]] <- F
    #for(j in 1:length(indexBands)){
    #  input[[bandNames[j]]][input[[qaBands[i]]] %in% look[qa == "bad", v]] <- NA #screen out qa values that are invalid (1)
    #}
  }
  
  # aerosol filter
  
  
  return(data.frame(valid,calcIndex(input, s=sat, b=bandNames)))
}


extract_cps <- function(bocpd_mod,truecp,cpthresh,cpmin){
  # first get cps
  cps <- unique(c(bocpd_mod$cpInds))[-1]
  # eliminate ones in the training period
  if(any(cps<cpmin)){
    cps <- cps[cps>=cpmin]
  }
  if(any(cps>(truecp+5))){
    cps <- cps[cps<=(truecp+5)]
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
#   if(any(cps>(truecp+5))){
#     cps <- cps[cps<=(truecp+5)]
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

# calculate_priors <- function(ptype,response_list,
#                              hls_sorted,change_date_p,
#                              npix,nsample=10000){
#   # select a sample of data to get prior info
#   pixels_sample <- sample(1:npix,nsample)
#   
#   sample_data_list <- lapply(pixels_sample,FUN=function(n){
#     datapts <- find_scale(pixel_num=n,response_list,
#                           hls_sorted,change_date_p,ptype)
#     return(datapts)
#   })
#   
#   # # calculate priors based on training data
#   # if(ptype==1){
#   #   get_priors <- get_priors1
#   # }
#   # if(ptype==2){
#   #   get_priors <- get_priors2
#   # }
#   # if(ptype==3){
#   #   get_priors <- get_priors3
#   # }
#   priors <- get_priors(data_list=sample_data_list,d=length(response_list),k=4)
#   return(priors)
# }
# 
# d=2
# k=4
# B0 <- matrix(c(1,0.1,0.2,0,0.5,-0.1,-0.2,0),4,2)
# Lambda0 <- diag(10,k)
# nu0 <- 50
# V0 <- (nu0-d-1)*matrix(c(1,0.7,0.7,1),2,2)


# sample parameters and data using true hyperparameters
# generate n time series with nt observed time points in each
# n <- 100
# nt <- 5*30
# tfull <- seq(0,5*365,length.out=nt)
# X <- cbind(1,sin(2*pi*tfull/365),cos(2*pi*tfull/365),tfull/(max(tfull)))
# 
# # sample parameters
# Sigmas <- replicate(n,riwish(nu0,V0),simplify=F)
# Betas <- lapply(Vs,function(x){
#   rmatnorm(1,B0,solve(Lambda0),x)
# })
# 
# # sample data
# Ys <- mapply(function(b,v){
#   X%*%b+mvrnorm(nt,rep(0,2),v)
# },b=Betas,v=Sigmas,SIMPLIFY = F)
# 
# # put data into format that the function expects
# Ys1 <- t(sapply(Ys,function(x){
#   c(x[,1])
# }))
# 
# Ys2 <- t(sapply(Ys,function(x){
#   c(x[,2])
# }))
# 
# response_list <- list(Ys1,Ys2)
# 
# 
# ####### Get hyperparameter estimates from the fixed, known parameters Beta and Sigma
# # calculate necessary quantities
# logdetV <- sum(sapply(Sigmas,function(x){
#   log(det(x))
# }),na.rm=T)
# meanV <- matrix(rowMeans(sapply(Sigmas,function(x){
#   c(x)
# }),na.rm=T),d,d)
# # function for the d/dnu loglikelihood
# f_nu <- function(nu,d,n,logdetV,meanV){
#   val1 <- (n/2)*log(det((nu-d-1)*meanV))
#   val2 <- -n*d*log(2)/2
#   val3 <- -0.5*logdetV
#   val4 <- -(n/2)*sum(sapply(1:d,function(j){
#     digamma(nu/2+(1-j)/2)
#   }))
#   val1+val2+val3+val4
# }
# # find nuhat by getting the root of f_nu
# root_result <- uniroot(f_nu,interval=c(d+1.01,200),d=d,n=n,logdetV=logdetV,meanV=meanV)
# 
# nuhat <- root_result$root
# 
# # estimate V
# Vhat <- meanV*(nuhat-d-1)
# # estimate B
# B <- t(sapply(Betas,function(x){c(x)}))
# Bhat <- matrix(colMeans(B,na.rm=T),k,d)
# 
# # estimate Lambdainv and Lambda
# LambdaInvhat <- (1/d)*matrix(colMeans(t(mapply(function(Bi,Sigmai,Bhat){
#   val1 <- Bi-Bhat
#   val1%*%solve(Sigmai)%*%t(val1)
# },Bi=Betas,Sigmai=Sigmas,MoreArgs = list(Bhat=Bhat))),na.rm=T),k,k)
# 
# Lambdahat <- solve(LambdaInvhat)
# 
# # compare estimates to true:
# Bhat
# B0
# Vhat
# V0
# Lambdahat
# Lambda0
# nuhat
# nu0
# 
# # Now get hyperparameter estimates from the data
# # Fit multivariate regression via gls to estimate Betahat and Sigmahat for each time series
# # Then estimate hyperparameters from Betahat and Sigmahat
# # Outputs the final hyperparameter estimates in "priors"
# # this takes a little bit 
# priors <- calculate_priors(3,response_list,hls_sorted=tfull,change_date_p=max(tfull),npix=n,nsample=n)
# 
# # compare estimates to true
# priors$B
# B0
# priors$V
# V0
# priors$Lambda
# Lambda0
# priors$nu
# nu0
