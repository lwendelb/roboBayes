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

get_bocpdms_metrics <- function(data,mincp=50){
  lapply(data,function(datai){
    lapply(datai,function(x){
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
        metrics <- get_metrics_from_cps(cps,cpd,tol=5)
        return(c(metrics,list(tm=t1)))
      }else{return(list(TP=NA,FP=NA,Precision=NA,Recall=NA,Fb=NA,latency=NA,tm=NA))}
      
    })
  })
  
}

get_other_metrics <- function(data_all,atype){
  lapply(data_all,function(datai){
    thing <- lapply(datai,function(cp_info){
      if(any(names(cp_info)=="bocpd_cps")){
        if(atype=="ccdc"){
          metrics <- get_metrics_from_cps(cp_info$ccdc_cps,cp_info$ccdc_cps+2,tol=5)
          t1 <- cp_info$ccdc_time
        }
        if(atype=="roboBayes"){
          metrics <- get_metrics_from_cps(cp_info$bocpd_cps,cp_info$bocpd_det,tol=5)
          t1 <- cp_info$bocpd_time
          
        }
        if(atype=="bocpd"){
          metrics <- get_metrics_from_cps(cp_info$bocpd2_cps,cp_info$bocpd2_det,tol=5)
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
    })
    return(thing)
  })
}


get_summary <- function(metrics,getsd=F){
  if(!getsd){
    TP <- sapply(metrics,function(metricsi){
      mean(sapply(metricsi,function(i){
        i$TP
      }),na.rm=T)
    })
    FP <- sapply(metrics,function(metricsi){
      mean(sapply(metricsi,function(i){
        i$FP
      }),na.rm=T)
    })
    Fb <- sapply(metrics,function(metricsi){
      mean(sapply(metricsi,function(i){
        i$Fb
      }),na.rm=T)
    })
    latency <- sapply(metrics,function(metricsi){
      mean(sapply(metricsi,function(i){
        i$latency
      }),na.rm=T)
    })
    tm <- 1000*sapply(metrics,function(metricsi){
      mean(sapply(metricsi,function(i){
        i$tm
      }),na.rm=T)/270
    })
  }else{
    TP <- sapply(metrics,function(metricsi){
      val <- sapply(metricsi,function(i){
        i$TP
      })
      sd(val,na.rm=T)/length(val)
    })
    FP <- sapply(metrics,function(metricsi){
      val <- sapply(metricsi,function(i){
        i$FP
      })
      sd(val,na.rm=T)/length(val)
    })
    Fb <- sapply(metrics,function(metricsi){
      val <- sapply(metricsi,function(i){
        i$Fb
      })
      sd(val,na.rm=T)/length(val)
    })
    latency <- sapply(metrics,function(metricsi){
      val <- sapply(metricsi,function(i){
        i$latency
      })
      sd(val,na.rm=T)/length(val)
    })
    tm <- sapply(metrics,function(metricsi){
      val <- sapply(metricsi,function(i){
        i$tm
      })
      sd(1000*val/270,na.rm=T)/length(val)
    })
  }
  
  return(cbind(TP,FP,Fb,latency,tm))
}


other_data <- list()
for(k in seq(1,11)){
  fn <- load(paste("results/ss",k,"sensitivity_v3.RData",sep=""))
  other_data <- c(other_data,list(all_results))
}

roboBayes_metrics <- get_other_metrics(other_data,atype="roboBayes")
roboBayes_mean <- get_summary(roboBayes_metrics,getsd=F)
roboBayes_sd <- get_summary(roboBayes_metrics,getsd=T)

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



library(xtable)


all_mns <- data.frame(roboBayes_mean)
all_se <- data.frame(roboBayes_sd)
colnames(all_mns) <- c("TP","FP","F-score","Latency","Time (ms)")
colnames(all_se) <- c("TP","FP","F-score","Latency","Time (ms)")

inds<- expand.grid(1:nrow(all_mns),1:ncol(all_mns))

mn_table <- all_mns
sd_table <- all_se
tab <- cbind(sim_settings,matrix(c(sapply(1:nrow(inds),function(y){
  paste(format(as.numeric(mn_table[inds[y,1],inds[y,2]]),nsmall=2,digits=1),
        "\\textsubscript{",
        format(as.numeric(sd_table[inds[y,1],inds[y,2]]),nsmall=2,digits=1),"}",sep="")
})),nrow(mn_table),ncol(mn_table)))

colnames(tab) <- c("$\\mu$","$\\Lambda_0$ scale","$\\nu_0$","$\\alpha$","$p_o$","L","TP","FP","F-score","Latency","Time (ms)")

cap <- "Average metrics (standard errors are in subscripts) for the simulation sensitivity
study applied to scenario 7.\\\\"
label <- "tab:sim_sensitivity_results"
print(xtable(tab, type = "latex",caption=cap,label=label),
      sanitize.text.function = identity,
      include.rownames = FALSE,
      caption.placement = "top",
      table.placement="H",
      hline.after=c(-1,0), 
      file = paste("Plots/ss_sensitivity_metrics_table_v3.tex",sep=""),
)


