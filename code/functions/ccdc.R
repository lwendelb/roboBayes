fit_stable <- function(y,x,mod=list()){
  XTXi <- mod$XTXi
  XTY <- mod$XTY
  nu <- mod$nu
  
  if(is.null(XTXi)){
    XTXi <- solve(t(x)%*%x)
    XTY <- t(x)%*%y
    nu <- nrow(x)
    return(list(XTXi=XTXi,XTY=XTY,nu=nu))
  }
  # update sufficient statistics
  XTXi_new <- XTXi-(XTXi%*%t(x)%*%(x)%*%XTXi)/c((1+(x)%*%XTXi%*%t(x)))
  XTY_new <- XTY + t(x)%*%y
  nu_new <- nu+1
  
  return(list(XTXi=XTXi_new,XTY=XTY_new,nu=nu_new))
}

runCCDC <- function(y,t,X,ts){
  X <- data.frame(X)
  ys <- lapply(y,function(x){x[1:ts]})
  ### fit initial stable period
  mods <- lapply(ys,fit_stable,x=as.matrix(X[1:ts,]))
  
  n <- length(y[[1]])
  k <- ncol(X)
  cps <- c(1)
  cpi <- 1
  outliers <- NULL
  ecount <- 0
  mod_storage <- list(mods)
  i=ts+1
  while(i <= n){
    sig <- mapply(function(mod,yy){
      ind_reset <- mod$nu + length(outliers)-(i-outliers)+1
      preds <- as.matrix(X[((i-mod$nu-length(outliers)):(i)),])%*%mod$XTXi%*%mod$XTY
      preds[ind_reset] <- yy[outliers]
      sqrt(sum((preds-yy[(i-mod$nu-length(outliers)):(i)])^2)/mod$nu)
    },mod=mods,yy=y)
    yi <- lapply(y,function(x){x[i]})
    #xi <- data.frame(x1=1, x2=sin(2*pi*t[i]/365), x3=cos(2*pi*t[i]/365), x4= (t[i]-min(t))/(max(t)-min(t)))
    xi <- as.matrix(X[i,])
    ypi <- lapply(mods,function(mod){
      xi%*%mod$XTXi%*%mod$XTY})
    ei <- mapply(function(yi_,ypi_,sig_){
      (yi_-ypi_)/sig_
    },yi_=yi,ypi_=ypi,sig_=sig)
    
    # check for change points
    cond1 <- (mean(abs(ei))/3) < 1
    if(cond1){
      # if a point is within bounds, reduce ecount and add to model
      ecount <- 0
      
      mods <- mapply(function(yy,md,xi){
        fit_stable(y=yy,xi,md)
      },yy=yi,md=mods,MoreArgs=list(xi=xi),SIMPLIFY=F )
        
    } else{ # if a point is not within bounds, increment ecount
      ecount <- ecount + 1
      outliers <- c(outliers,i)
    }
    if(ecount==3){
      # record cp
      cp <- t[i-2]
      cpi <- i-2
      cps <- c(cps,cpi)
      
      # reinitialize
      if((length(t)-i-2) >= (ts)){
        mods <- lapply(y, function(yy){
          fit_stable(y=yy[(i-2):(i+ts-3)],x=as.matrix(X[(i-2):(i+ts-3),]))
        } )
      }
      else{
        i <- n
      }
      
      # increment i
      i <- i+ts+1
      ecount=0
      outliers <- NULL
    }else{
      i <- i + 1
    }
    #mod_storage <- c(mod_storage,list(mods))
    
  }
  return(list(cps=cps))
}

# results_ccdc <- runCCDC(Y,t,mods)
# results_ccdc$mod_storage[[1]]
# results_ccdc$cps
# length(results_ccdc$mod_storage)
# plot(t,Y[[1]],ylim=c(0,10))
# abline(v=results_ccdc$cps)

runBFAST <- function(y,t,tfreq=20,start_cp=c(3,1)){
  y <- as.matrix(y)
  d <- ncol(y)
  cps <- c()
  cp_det <- c()
  for(i in ((start_cp[1]-1)*tfreq+start_cp[2]):nrow(y)){
    ii <- i-tfreq
    start_monitor <- c(floor(i/tfreq)+1,2*(ii%%tfreq))
    #print(start_monitor)
    ti <- t[1:i]
    y.ts <- ts(y[1:i,],deltat=t[2]-t[1],frequency=tfreq)
    xi <- cbind(x1=1, x2=sin(2*pi*ti/365), x3=cos(2*pi*ti/365), x4= (ti)/1000)
    
    y2 <- c(t(y.ts))
    y2.ts <- ts(y2,deltat=t[2]-t[1],frequency=d*tfreq)
    x2i1 <- cbind(xi,0,0,0,0)
    x2i2 <- cbind(0,0,0,0,xi)
    X2 <- matrix(c(t(cbind(x2i1,x2i2))),d*length(ti),8,byrow=T)
    xreg <- X2
    bf <- bfastmonitor(data=y2.ts,start=c(6,1),history=c(1,1),formula = response ~ xreg - 1,order=0)
    
    cp <- floor((bf$breakpoint-1)*tfreq+1)
    
    #print(c(cp,i))
    if(!is.na(cp)){
      if(length(cps>0)){
        if(!any(cps==cp)){
          cps <- c(cps,cp)
          cp_det <- c(cp_det,i)
          start_monitor <- cp+1
        }
      }else{
          cps <- c(cps,cp)
          cp_det <- c(cp_det,i)
          start_monitor <- cp+1
      }
    }
  }
 return(list(cps=cps,cp_det=cp_det))
}


# fit_stable <- function(y,t,tchange){
#   
#   # pull out stable observations
#   ind_stable <- which(t <= tchange)
#   y_stable <- y[ind_stable]
#   t_stable <- t[ind_stable]
#   
#   # construct stable 
#   # fit TS model
#   x1 <- 1
#   x2 <- sin(2*pi*t_stable/365)
#   x3 <- cos(2*pi*t_stable/365)
#   x4 <- (t_stable-min(t))/(max(t)-min(t))
#   X_stable <- data.frame(x1,x2,x3,x4)
#   
#   lm_mod <- lm(y_stable ~.-1, X_stable,na.action=na.omit)
#   
#   return(lm_mod)
# }
# 
# runCCDC <- function(y,t,X){
#   
#   ### fit initial stable period
#   mods <- lapply(y,fit_stable,t=t,tchange=t[12])
#   
#   n <- length(y[[1]])
#   cps <- c()
#   cpi <- 1
#   ecount <- 0
#   mod_storage <- list(mods)
#   i=13
#   while(i <= n){
#     sig <- lapply(mods,sigma)
#     yi <- lapply(y,function(x){x[i]})
#     xi <- data.frame(x1=1, x2=sin(2*pi*t[i]/365), x3=cos(2*pi*t[i]/365), x4= (t[i]-min(t))/(max(t)-min(t)))
#     ypi <- lapply(mods,function(mod){
#       predict(mod,xi)})
#     ei <- mapply(function(yi_,ypi_,sig_){
#       (yi_-ypi_)/sig_
#     },yi_=yi,ypi_=ypi,sig_=sig)
#     
#     # check for change points
#     cond1 <- (mean(abs(ei))/3) < 1
#     if(cond1){
#       # if a point is within bounds, reduce ecount and add to model
#       ecount <- 0
#       
#       mods <- lapply(y, function(yy){
#         fit_stable(y=yy[(cpi):(i)],t=t[cpi:i],tchange=t[i])
#       } )
#       
#     } else{ # if a point is not within bounds, increment ecount
#       ecount <- ecount + 1
#     }
#     if(ecount==3){
#       # record cp
#       cp <- t[i-2]
#       cpi <- i-2
#       cps <- c(cps,cpi)
#       
#       # reinitialize
#       if((length(t)-i-2) >= 15){
#         mods <- lapply(y, function(yy){
#           fit_stable(y=yy[(i-2):(i+9)],t=t[(i-2):(i+9)],tchange=t[(i+9)])
#         } )
#       }
#       else{
#         i <- n
#       }
#       
#       # increment i
#       i <- i+12
#       ecount=0
#     }else{
#       i <- i + 1
#     }
#     mod_storage <- c(mod_storage,list(mods))
#     
#   }
#   return(list(cps=cps,mod_storage=mod_storage))
# }
# 
# # results_ccdc <- runCCDC(Y,t,mods)
# # results_ccdc$mod_storage[[1]]
# # results_ccdc$cps
# # length(results_ccdc$mod_storage)
# # plot(t,Y[[1]],ylim=c(0,10))
# # abline(v=results_ccdc$cps)
# 
# runBFAST <- function(y,t,tfreq=20,start_cp=c(3,1)){
#   y <- as.matrix(y)
#   d <- ncol(y)
#   cps <- c()
#   cp_det <- c()
#   for(i in ((start_cp[1]-1)*tfreq+start_cp[2]):nrow(y)){
#     ii <- i-tfreq
#     start_monitor <- c(floor(i/tfreq)+1,2*(ii%%tfreq))
#     #print(start_monitor)
#     ti <- t[1:i]
#     y.ts <- ts(y[1:i,],deltat=t[2]-t[1],frequency=tfreq)
#     xi <- cbind(x1=1, x2=sin(2*pi*ti/365), x3=cos(2*pi*ti/365), x4= (ti)/1000)
#     
#     y2 <- c(t(y.ts))
#     y2.ts <- ts(y2,deltat=t[2]-t[1],frequency=d*tfreq)
#     x2i1 <- cbind(xi,0,0,0,0)
#     x2i2 <- cbind(0,0,0,0,xi)
#     X2 <- matrix(c(t(cbind(x2i1,x2i2))),d*length(ti),8,byrow=T)
#     xreg <- X2
#     bf <- bfastmonitor(data=y2.ts,start=c(6,1),history=c(1,1),formula = response ~ xreg - 1,order=0)
#     
#     cp <- floor((bf$breakpoint-1)*tfreq+1)
#     
#     #print(c(cp,i))
#     if(!is.na(cp)){
#       if(length(cps>0)){
#         if(!any(cps==cp)){
#           cps <- c(cps,cp)
#           cp_det <- c(cp_det,i)
#           start_monitor <- cp+1
#         }
#       }else{
#         cps <- c(cps,cp)
#         cp_det <- c(cp_det,i)
#         start_monitor <- cp+1
#       }
#     }
#   }
#   return(list(cps=cps,cp_det=cp_det))
# }


# pt <- df[df$pointid=="g8p14_2",]
# #pt$t <- pt$t %% 365
# library(GPfit)
# cor_fun <- function(d){
#   l <- 2
#   exp((-2*sin(d)^2)/l^2)
# }
# 
# Yj <- rnorm(length(pt$t),pt$t,0.01)
# Yj <- (Yj-min(Yj))/(max(Yj)-min(Yj))
# X <- cbind(cos(2*pi*pt$t/365)/2+0.5,Yj)[1:90,]
# #X <- apply(1,)
# Y <- (pt$index-min(pt$index[1:90]))/(max(pt$index[1:90])-min(pt$index[1:90]))
# gpmod <- GP_fit(X=X[,1:2],Y=pt$index[1:90],corr = list(type = "exponential",power=1.9),nug_thres=10)
# tnew <- min(pt$t):max(pt$t)
# Xnew <- cbind(cos(2*pi*tnew/365)/2+0.5,(tnew-min(tnew))/(max(tnew)-min(tnew)))
# preds <- predict(gpmod,Xnew[,1:2])
# #plot.GP(gpmod)
# plot(pt$t,pt$index)
# lines(tnew,preds$Y_hat)
# lines(tnew,preds$Y_hat+2*sqrt(preds$MSE))
# 
# #predict(gpmod)
# 
# 
# library(mlegp)
# Y <- pt$index[1:90]
# 
# X <- cbind(cos(2*pi*pt$t/365)/2+0.5,pt$t/20000)[1:90,]
# Xnew <- cbind(cos(2*pi*tnew/365)/2+0.5,tnew/20000)
# 
# 
# gpmod <- mlegp(X=X[,1],Z=Y,nugget=.002,nugget.known=T,parallel=F)
# gpmod$beta
# preds <- predict.gp(gpmod,newData=as.matrix(Xnew[,1]),se.fit=T)
# plot(tnew,preds$fit,type="l",ylim=c(0,1))
# points(pt$t,pt$index)
# lines(tnew,preds$fit+2*preds$se.fit+2*sqrt(.002))
# lines(tnew,preds$fit-2*preds$se.fit-2*sqrt(.002))
# 
# 
# library(GauPro)
# X <- cbind(1,sin(2*pi*pt$t/365),cos(2*pi*pt$t/365))[1:90,]
# 
# tnew <- min(pt$t):max(pt$t)
# Xnew <- cbind(1,sin(2*pi*tnew/365)/2+0.5,cos(2*pi*tnew/365)/2+0.5)
# trend <- trend_LM$new(D=3)
# gp <- GauPro(X[,3],Y,parallel=F)
# preds <- predict(gp,XX=Xnew,se.fit=T)
# #plot(tnew,preds$mean)
# #plot(tnew,preds$se)
# plot(tnew,preds$mean,type='l',ylim=c(0,1))
# lines(tnew,preds$mean+2*preds$se)
# lines(tnew,preds$mean-2*preds$se)
# points(pt$t,pt$index)
# 
# library(laGP)
# laGP(Xref=)
# 
# X <- cbind(1,sin(2*pi*pt$t/365),cos(2*pi*pt$t/365))[1:90,]
# 
# lamod <- likfit(coords = cbind(rnorm(90,pt$t[1:90],.1),X[,2]), data = as.matrix(Y),
#        trend = ~ X[,2:3],ini=c(var(Y),1),nugget=.01,fix.kappa=F,fix.lambda=F)
# krige.conv(lamod)
# krige.conv(parana, loc = parana.gr, krige = KC, output = OC)
# predict(lamod)
# predGP(lamod,XX=)
