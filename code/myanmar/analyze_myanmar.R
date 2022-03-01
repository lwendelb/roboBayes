get_metrics_from_cps <- function(cps,cpd,true_cp=181,L=5,p=0,n=270,tol=2,beta=1){
  true_cp_tol <- sort(c(true_cp-(min(1,tol):tol),true_cp,true_cp+(min(1,tol):tol)))
  true_cp_tol <- true_cp_tol[(true_cp_tol<=(n))&(true_cp_tol>(L+p))]
  
  TP_count  <- 0*cps
  FP_count <- 0*cps
  
  latency_TP <- 0*cps
  if(length(cps)>=1){
    for(i in 1:length(cps)){
      TP_count[i] <- any(cps[i]==true_cp_tol)
      FP_count[i] <- all(cps[i]!=true_cp_tol)
    }
    
    if(any(TP_count>=1)){
      latency_TP <- min((cpd-cps)[TP_count>=1])
      
    }else{
      latency_TP <- NA
    }
    
    TP <- as.numeric(sum(TP_count)>0)
    FP <- sum(FP_count)
  }else{
    TP <- 0
    FP <- 0
    latency_TP <- NA
  }
  
  Precision <- (TP+1)/(TP+FP+1)
  Recall <- (TP+1)/(length(true_cp)+1)
  
  Fb <- (1+beta^2)*Precision*Recall/((beta^2)*Precision + Recall)
  
  Ntot <- n-tol*2-1
  
  return(list(TP=TP,FP=FP,Precision=Precision,Recall=Recall,Fb=Fb,latency_TP=latency_TP))
  
}

get_bocpdms_metrics <- function(data,truecps,mincp=12){
  mapply(function(x,truecp){
    if(any(names(x)=="rbocpdms" | names(x)=="bocpdms")){
      if(any(names(x)=="rbocpdms")){
        x1 <- x$rbocpdms
        t1 <- x$rbocpdms_time
      }
      if(any(names(x)=="bocpdms")){
        x1 <- x$bocpdms
        t1 <- x$bocpdms_time
      }
      # first get cps
      cps <- unique(unlist(x1))
      # eliminate ones in the training period
      if(any(cps<mincp)){
        cps <- cps[cps>=mincp]
      }
      if(any(cps>(truecp+5))){
        cps <- cps[cps<=(truecp+5)]
      }
      # get latency
      cpd <- 0*cps
      if(length(cps)>0){
        for(cpi in 1:length(cps)){
          cpd[cpi] <- which(sapply(x1,function(xi){
            any(cps[cpi]==xi)
          }))[1]+1
        }
      }
      
      # now get TP and FP metrics
      metrics <- get_metrics_from_cps(cps,cpd,true_cp=truecp,tol=5)
      return(c(metrics,list(tm=t1)))
    }else{return(list(TP=NA,FP=NA,Precision=NA,Recall=NA,Fb=NA,latency=NA,tm=NA))}
    
  },x=data,truecp=truecps,SIMPLIFY = F)
  
}

get_other_metrics <- function(datai,truecps,atype){
    thing <- mapply(function(cp_info,truecp){
      if(any(names(cp_info)=="ccdc_cps")){
        if(atype=="ccdc"){
          metrics <- get_metrics_from_cps(cp_info$ccdc_cps,cp_info$ccdc_cps+2,true_cp=truecp,tol=5)
          t1 <- cp_info$ccdc_time
        }
        if(atype=="bocpdod"){
          metrics <- get_metrics_from_cps(cp_info$bocpd_cps,cp_info$bocpd_det,true_cp=truecp,tol=5)
          t1 <- cp_info$bocpd_time
          
        }
        if(atype=="bocpd"){
          metrics <- get_metrics_from_cps(cp_info$bocpd2_cps,cp_info$bocpd2_det,true_cp=truecp,tol=5)
          t1 <- cp_info$bocpd2_time
          
          
        }
        #bfast_metrics <- get_metrics_from_cps(cp_info$bfast_cps,cp_info$bfast_detected,tol=2)
        
        
        
        # bocpd_metrics <- mapply(function(x,y){
        #   get_metrics_from_cps(x,y,tol=2)
        # },cp_info$bocpd_cps,cp_info$bocpd_det)
        #bfast_latency <- cp_info$bfast_detected-cp_info$bfast_cps
        
        # bocpd_latency <- mapply(function(x1,x2){
        #   ifelse(length(x2-x1)>=1,x2-x1,NA)
        # },cp_info$bocpd_cps,cp_info$bocpd_det)
        # 
        # bocpd2_latency <- mapply(function(x1,x2){
        #   ifelse(length(x2-x1)>=1,x2-x1,NA)
        # },cp_info$bocpd2_cps,cp_info$bocpd2_det)
        
        #bocpd_metrics <- get_metrics_from_cps(cp_info$bocpd_cps)
        #bocpd_latency <- cp_info$bocpd_det - cp_info$bocpd_cps
        return(c(metrics,list(tm=t1)))
        # return(list(ccdc_TP=ccdc_metrics$TP,ccdc_FP=ccdc_metrics$FP,
        #             bfast_TP=bfast_metrics$TP,bfast_FP=bfast_metrics$FP,
        #             bocpd_metrics=bocpd_metrics,
        #             bocpd2_metrics=bocpd2_metrics,
        #             bocpd_TP=bocpd_metrics$TP,bocpd_FP=bocpd_metrics$FP,
        #             bocpd2_TP=bocpd2_metrics$TP,bocpd2_FP=bocpd2_metrics$FP,
        #             bfast_latency=bfast_latency,
        #             bocpd_latency=bocpd_latency,
        #             bocpd2_latency=bocpd2_latency,
        #             bfast_latency_TP=bfast_metrics$latency_TP,
        #             bocpd_latency_TP=bocpd_metrics$latency_TP,
        #             bocpd2_latency_TP=bocpd2_metrics$latency_TP))
      }
      else{return(list(TP=NA,FP=NA,Precision=NA,Recall=NA,Fb=NA,latency=NA,tm=NA))}
    },cp_info=datai,truecp=truecps,SIMPLIFY = F)
    return(thing)
}


get_summary <- function(metricsi,getsd=F){
  if(!getsd){
    TP <-  mean(sapply(metricsi,function(i){
        i$TP
      }),na.rm=T)
    
    FP <- mean(sapply(metricsi,function(i){
        i$FP
      }),na.rm=T)
    
    Fb <- mean(sapply(metricsi,function(i){
        i$Fb
      }),na.rm=T)
    
    latency <-mean(sapply(metricsi,function(i){
        i$latency
      }),na.rm=T)
    
    tm <- mean(length(dfl)*1000*sapply(metricsi,function(i){
        i$tm
      })/sum(sapply(dfl,function(x){nrow(x)})),na.rm=T)
    
  }else{
    TP <- sd(sapply(metricsi,function(i){
        i$TP
      }),na.rm=T)
    
    FP <- sd(sapply(metricsi,function(i){
        i$FP
      }),na.rm=T)
    
    Fb <- sd(sapply(metricsi,function(i){
        i$Fb
      }),na.rm=T)
    
    latency <- sd(sapply(metricsi,function(i){
        i$latency
      }),na.rm=T)
    
    tm <- sd(length(dfl)*1000*sapply(metricsi,function(i){
        i$tm
      })/sum(sapply(dfl,function(x){nrow(x)})),na.rm=T)
    
  }
  
  return(cbind(TP,FP,Fb,latency,tm))
}

load("Myanmar/myanmar_rbocpdms_p5_p5_100.RData")
data_rbocpdms <- all_results
load("Myanmar/myanmar_bocpdms_1_1_100.RData")
data_bocpdms <- all_results


load("Myanmar/myanmar_v2_rbocpdms.RData")
data_rbocpdms <- all_results
load("Myanmar/myanmar_v2_bocpdms.RData")
data_bocpdms <- all_results
load("v2_results/myanmar_v2.RData")
data_other <- all_results

load("v3_results/myanmar_v3_rbocpdms.RData")
data_rbocpdms <- all_results
load("v3_results/myanmar_v3_bocpdms.RData")
data_bocpdms <- all_results
load("v3_results/myanmar_v3.RData")
data_other <- all_results

truecps <- lapply(dfl,function(dfi){
  dfi$ti[which(dfi$date_dist<=dfi$date)[1]]
})

dist_perc <- sapply(dfl,function(dfi){
  dfi$dist_perc[1]
})

metrics_bocpdms <- get_bocpdms_metrics(data_bocpdms,truecps)
bocpdms_mean <- get_summary(metrics_bocpdms)

metrics_rbocpdms <- get_bocpdms_metrics(data_rbocpdms,truecps)
rbocpdms_mean <- get_summary(metrics_rbocpdms)

metrics_ccd <- get_other_metrics(data_other,truecps,atype="ccdc")
ccd_mean <- get_summary(metrics_ccd)
metrics_bocpd <- get_other_metrics(data_other,truecps,atype="bocpd")
bocpd_mean <- get_summary(metrics_bocpd)
metrics_bocpdod <- get_other_metrics(data_other,truecps,atype="bocpdod")
bocpdod_mean <- get_summary(metrics_bocpdod)

mn_table <- rbind(ccd_mean,bocpdms_mean,rbocpdms_mean,bocpd_mean,bocpdod_mean)

mn_table <- apply(mn_table,c(1,2),function(x){
  format(x,digits=1,nsmall=2)
})

tab <- rbind(c("TP","FP","F-score","Latency","Time (ms)"),mn_table)

rownames(tab) <- c("Metric","CCD","BOCPDMS","rBOCPDMS","BOCPD","RobOBayes")
#colnames(tab) <- c("Metric","CCD","BOCPDMS","rBOCPDMS","BOCPD","BOCPD-OD")

library(xtable)
cap <- c("Summary statistics for detection of the first change point in annotated Myanmar deforestation data.")
label <- "tab:myanmar_results"
print(xtable(tab, type = "latex",caption=cap,label=label),
      sanitize.text.function = identity,
      include.colnames = FALSE,
      caption.placement = "top",
      table.placement="H",
      hline.after=c(-1,1), file = paste("myanmar_results_table_v3.tex",sep=""),
)




######### subsets based on change size
load("v3_results/myanmar_v3_rbocpdms.RData")
data_rbocpdms <- all_results
load("v3_results/myanmar_v3_bocpdms.RData")
data_bocpdms <- all_results
load("v3_results/myanmar_v3.RData")
data_other <- all_results


ind_low <- dist_perc < 0.5
ind_med <- dist_perc >=0.5 & dist_perc < 0.9
ind_high <- dist_perc >= 0.9

data_other_low <- data_other[ind_low]
data_other_med <- data_other[ind_med]
data_other_high <- data_other[ind_high]

true_cp_list <- list(truecps[ind_low],truecps[ind_med],truecps[ind_high])
data_other_list <- list(data_other_low,data_other_med,data_other_high)
data_bocpdms_list <- list(data_bocpdms[ind_low],data_bocpdms[ind_med],data_bocpdms[ind_high])
data_rbocpdms_list <- list(data_rbocpdms[ind_low],data_rbocpdms[ind_med],data_rbocpdms[ind_high])

metrics_ccdc <- mapply(function(dat,tcp){
  get_other_metrics(dat,tcp,atype="ccdc")
},dat=data_other_list,tcp=true_cp_list)
ccd_mean <- do.call(rbind,lapply(metrics_ccdc,function(x){
  get_summary(x)
}))

metrics_bocpd <- mapply(function(dat,tcp){
  get_other_metrics(dat,tcp,atype="bocpd")
},dat=data_other_list,tcp=true_cp_list)
bocpd_mean <- do.call(rbind,lapply(metrics_bocpd,function(x){
  get_summary(x)
}))

metrics_bocpdod <- mapply(function(dat,tcp){
  get_other_metrics(dat,tcp,atype="bocpdod")
},dat=data_other_list,tcp=true_cp_list)
bocpdod_mean <- do.call(rbind,lapply(metrics_bocpdod,function(x){
  get_summary(x)
}))

metrics_bocpdms <- mapply(function(dat,tcp){
  get_bocpdms_metrics(dat,tcp)
},dat=data_bocpdms_list,tcp=true_cp_list)
bocpdms_mean <- do.call(rbind,lapply(metrics_bocpdms,function(x){
  get_summary(x)
}))

metrics_rbocpdms <- mapply(function(dat,tcp){
  get_bocpdms_metrics(dat,tcp)
},dat=data_rbocpdms_list,tcp=true_cp_list)
rbocpdms_mean <- do.call(rbind,lapply(metrics_rbocpdms,function(x){
  get_summary(x)
}))



mn_table <- rbind(ccd_mean,bocpdms_mean,
                  rbocpdms_mean,
                  bocpd_mean,bocpdod_mean)

mn_table <- apply(mn_table,c(1,2),function(x){
  format(x,digits=1,nsmall=2)
})

mn_table <- cbind(rep(c("small","med","large"),5),mn_table)

tab <- data.frame(mn_table)

tab <- tab %>% arrange(V1)

tab <- cbind(rep(c("CCD","BOCPDMS","rBOCPDMS","BOCPD","BOCPD-OD"),3),tab)
#colnames(tab) <- c("Metric","CCD","BOCPDMS","rBOCPDMS","BOCPD","BOCPD-OD")

colnames(tab) <- c("Method","Disturbance Size","TP","FP","F-score","Latency","Time (ms)")


library(xtable)
cap <- c("Summary statistics for detection of the first change point in annotated Myanmar deforestation data.")
label <- "tab:myanmar_results_v3"
print(xtable(tab, type = "latex",caption=cap,label=label),
      sanitize.text.function = identity,
      include.colnames = TRUE,
      include.rownames = FALSE,
      caption.placement = "top",
      table.placement="H",
      hline.after=c(-1,0,5,10), file = paste("myanmar_results_table_v3.tex",sep=""),
)


