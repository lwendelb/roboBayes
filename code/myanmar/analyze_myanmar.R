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
        if(atype=="roboBayes"){
          metrics <- get_metrics_from_cps(cp_info$bocpd_cps,cp_info$bocpd_det,true_cp=truecp,tol=5)
          t1 <- cp_info$bocpd_time
          
        }
        if(atype=="bocpd"){
          metrics <- get_metrics_from_cps(cp_info$bocpd2_cps,cp_info$bocpd2_det,true_cp=truecp,tol=5)
          t1 <- cp_info$bocpd2_time
          
          
        }
        
        return(c(metrics,list(tm=t1)))
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

load("results/myanmar_v3_rbocpdms.RData")
data_rbocpdms <- all_results
load("results/myanmar_v3_bocpdms.RData")
data_bocpdms <- all_results
load("results/myanmar_v3.RData")
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
metrics_roboBayes <- get_other_metrics(data_other,truecps,atype="roboBayes")
roboBayes_mean <- get_summary(metrics_roboBayes)

mn_table <- rbind(ccd_mean,bocpdms_mean,rbocpdms_mean,bocpd_mean,roboBayes_mean)

mn_table <- apply(mn_table,c(1,2),function(x){
  format(x,digits=1,nsmall=2)
})

tab <- rbind(c("TP","FP","F-score","Latency","Time (ms)"),mn_table)

rownames(tab) <- c("Metric","CCD","BOCPDMS","rBOCPDMS","BOCPD","roboBayes")

library(xtable)
cap <- c("Summary statistics for detection of the first change point in annotated Myanmar deforestation data.")
label <- "tab:myanmar_results"
print(xtable(tab, type = "latex",caption=cap,label=label),
      sanitize.text.function = identity,
      include.colnames = FALSE,
      caption.placement = "top",
      table.placement="H",
      hline.after=c(-1,1), file = paste("Plots/myanmar_results_table_v3.tex",sep=""),
)

