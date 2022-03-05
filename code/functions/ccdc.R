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

