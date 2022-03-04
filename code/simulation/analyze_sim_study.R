bocpdms_data <- list()
for(k in c(seq(1,4),9)){
  fn <- load(paste("results/ss",k,"bocpdms_default.RData",sep=""))
  bocpdms_data <- c(bocpdms_data,list(all_results))
}

rbocpdms_data <- list()
for(k in seq(1,9)){
  fn <- load(paste("results/ss",k,"rbocpdms_default.RData",sep=""))
  rbocpdms_data <- c(rbocpdms_data,list(all_results))
}
other_data <- list()
for(k in seq(1,9)){
  fn <- load(paste("results/ss",k,"sim_v3.RData",sep=""))
  other_data <- c(other_data,list(all_results))
}


bocpdms_metrics <- get_bocpdms_metrics(bocpdms_data)
bocpdms_mean <- get_summary(bocpdms_metrics,getsd=F)
bocpdms_sd <- get_summary(bocpdms_metrics,getsd=T)

rbocpdms_metrics <- get_bocpdms_metrics(rbocpdms_data)
rbocpdms_mean <- get_summary(rbocpdms_metrics,getsd=F)
rbocpdms_sd <- get_summary(rbocpdms_metrics,getsd=T)

ccdc_metrics <-   get_other_metrics(other_data,atype="ccdc")
ccdc_mean <- get_summary(ccdc_metrics,getsd=F)
ccdc_sd <- get_summary(ccdc_metrics,getsd=T)


bocpdod_metrics <- get_other_metrics(other_data,atype="bocpdod")
bocpdod_mean <- get_summary(bocpdod_metrics,getsd=F)
bocpdod_sd <- get_summary(bocpdod_metrics,getsd=T)


bocpd_metrics <-get_other_metrics(other_data,atype="bocpd")
bocpd_mean <- get_summary(bocpd_metrics,getsd=F)
bocpd_sd <- get_summary(bocpd_metrics,getsd=T)




library(xtable)
sim_settings <- cbind(c(c(0.1,0.2,rep(c(0.1,0.2),3))),
                      c(rep(c(0,0,0.9,0.9),2)),
                      c(rep(0,4),rep(0.05,4)))
#combine simulations 1-4
i <- 1
TP_mns <- (cbind(1:4,(ccdc_mean[1:4,i]),bocpdms_mean[1:4,i],
                    (rbocpdms_mean[1:4,i]),bocpd_mean[1:4,i],
                    bocpdod_mean[1:4,i]))
TP_sds <- (cbind(1:4,(ccdc_sd[1:4,i]),bocpdms_sd[1:4,i],
                   (rbocpdms_sd[1:4,i]),bocpd_sd[1:4,i],
                   bocpdod_sd[1:4,i]))
i <- 2
FP_mns <- (cbind(1:4,(ccdc_mean[1:4,i]),bocpdms_mean[1:4,i],
                 (rbocpdms_mean[1:4,i]),bocpd_mean[1:4,i],
                 bocpdod_mean[1:4,i]))
FP_sds <- (cbind(1:4,(ccdc_sd[1:4,i]),bocpdms_sd[1:4,i],
                 (rbocpdms_sd[1:4,i]),bocpd_sd[1:4,i],
                 bocpdod_sd[1:4,i]))
i <- 3
Fb_mns <- (cbind(1:4,(ccdc_mean[1:4,i]),bocpdms_mean[1:4,i],
                 (rbocpdms_mean[1:4,i]),bocpd_mean[1:4,i],
                 bocpdod_mean[1:4,i]))
Fb_sds <- (cbind(1:4,(ccdc_sd[1:4,i]),bocpdms_sd[1:4,i],
                 (rbocpdms_sd[1:4,i]),bocpd_sd[1:4,i],
                 bocpdod_sd[1:4,i]))
i <- 4
latency_mns <- (cbind(1:4,(ccdc_mean[1:4,i]),bocpdms_mean[1:4,i],
                 (rbocpdms_mean[1:4,i]),bocpd_mean[1:4,i],
                 bocpdod_mean[1:4,i]))
latency_sds <- (cbind(1:4,(ccdc_sd[1:4,i]),bocpdms_sd[1:4,i],
                 (rbocpdms_sd[1:4,i]),bocpd_sd[1:4,i],
                 bocpdod_sd[1:4,i]))

i <- 5
time_mns <- (cbind(1:4,(ccdc_mean[1:4,i]),bocpdms_mean[1:4,i],
                 (rbocpdms_mean[1:4,i]),bocpd_mean[1:4,i],
                 bocpdod_mean[1:4,i]))
time_sds <- (cbind(1:4,(ccdc_sd[1:4,i]),bocpdms_sd[1:4,i],
                 (rbocpdms_sd[1:4,i]),bocpd_sd[1:4,i],
                 bocpdod_sd[1:4,i]))


#combine simulations 5:8
i <- 1
TP_mns2 <- (cbind(5:8,(ccdc_mean[5:8,i]),bocpd_mean[5:8,i],
                 bocpdod_mean[5:8,i]))
TP_sds2 <- (cbind(5:8,(ccdc_sd[5:8,i]),bocpd_sd[5:8,i],
                 bocpdod_sd[5:8,i]))

TP_mns3 <- cbind(9,ccdc_mean[9,i],bocpdms_mean[5,i],
                 (rbocpdms_mean[9,i]),bocpd_mean[9,i],
                 bocpdod_mean[9,i])
TP_sds3 <- cbind(9,ccdc_sd[9,i],bocpdms_sd[5,i],
                 (rbocpdms_sd[9,i]),bocpd_sd[9,i],
                 bocpdod_sd[9,i])

i <- 2
FP_mns2 <- (cbind(5:8,(ccdc_mean[5:8,i]),bocpd_mean[5:8,i],
                 bocpdod_mean[5:8,i]))
FP_sds2 <- (cbind(5:8,(ccdc_sd[5:8,i]),bocpd_sd[5:8,i],
                 bocpdod_sd[5:8,i]))

FP_mns3 <- cbind(9,ccdc_mean[9,i],bocpdms_mean[5,i],
                 (rbocpdms_mean[9,i]),bocpd_mean[9,i],
                 bocpdod_mean[9,i])
FP_sds3 <- cbind(9,ccdc_sd[9,i],bocpdms_sd[5,i],
                 (rbocpdms_sd[9,i]),bocpd_sd[9,i],
                 bocpdod_sd[9,i])


i <- 3
Fb_mns2 <- (cbind(5:8,(ccdc_mean[5:8,i]),bocpd_mean[5:8,i],
                 bocpdod_mean[5:8,i]))
Fb_sds2 <- (cbind(5:8,(ccdc_sd[5:8,i]),bocpd_sd[5:8,i],
                 bocpdod_sd[5:8,i]))

Fb_mns3 <- cbind(9,ccdc_mean[9,i],bocpdms_mean[5,i],
                 (rbocpdms_mean[9,i]),bocpd_mean[9,i],
                 bocpdod_mean[9,i])
Fb_sds3 <- cbind(9,ccdc_sd[9,i],bocpdms_sd[5,i],
                 (rbocpdms_sd[9,i]),bocpd_sd[9,i],
                 bocpdod_sd[9,i])


i <- 4
latency_mns2 <- (cbind(5:8,(ccdc_mean[5:8,i]),bocpd_mean[5:8,i],
                      bocpdod_mean[5:8,i]))
latency_sds2 <- (cbind(5:8,(ccdc_sd[5:8,i]),bocpd_sd[5:8,i],
                      bocpdod_sd[5:8,i]))

latency_mns3 <- cbind(9,ccdc_mean[9,i],bocpdms_mean[5,i],
                      (rbocpdms_mean[9,i]),bocpd_mean[9,i],
                 bocpdod_mean[9,i])
latency_sds3 <- cbind(9,ccdc_sd[9,i],bocpdms_sd[5,i],
                      (rbocpdms_sd[9,i]),bocpd_sd[9,i],
                 bocpdod_sd[9,i])



i <- 5
time_mns2 <- (cbind(5:8,(ccdc_mean[5:8,i]),bocpd_mean[5:8,i],
                   bocpdod_mean[5:8,i]))
time_sds2 <- (cbind(5:8,(ccdc_sd[5:8,i]),bocpd_sd[5:8,i],
                   bocpdod_sd[5:8,i]))

time_mns3 <- cbind(9,ccdc_mean[9,i],bocpdms_mean[5,i],
                      (rbocpdms_mean[9,i]),bocpd_mean[9,i],
                      bocpdod_mean[9,i])
time_sds3 <- cbind(9,ccdc_sd[9,i],bocpdms_sd[5,i],
                      (rbocpdms_sd[9,i]),bocpd_sd[9,i],
                      bocpdod_sd[9,i])


#### set up simulation 9 stuff



make_table <- function(mn_table,sd_table,atype,scens=1:4){
  inds<- expand.grid(1:nrow(mn_table),2:ncol(mn_table))
  
  tab <- matrix(c(scens,sapply(1:nrow(inds),function(y){
    paste(format(mn_table[inds[y,1],inds[y,2]],nsmall=2,digits=1),
          "\\textsubscript{",
          format(sd_table[inds[y,1],inds[y,2]],nsmall=2,digits=1),"}",sep="")
  })),length(scens),ncol(mn_table))
  
  if(ncol(tab)==4){
    colnames(tab) <- c("Scenario","CCD","BOCPD","RobOBayes")
    
  }else{
    colnames(tab) <- c("Scenario","CCD","BOCPDMS","rBOCPDMS","BOCPD","RobOBayes")
    
  }

  print(xtable(tab, type = "latex",caption=atype),
        sanitize.text.function = identity,
        include.rownames = FALSE,
        hline.after=c(-1,0), file = paste("ss_",atype,scens[1],scens[length(scens)],"_table.tex",sep=""),
        )
}

methods <- c("CCD","BOCPDMS","rBOCPDMS","BOCPD","RobOBayes")
scens <- c(sapply(1:4,function(i){rep(i,5)}))
all_mns1 <- data.frame(scens,rep(methods,4),c(t(TP_mns[,-1])),c(t(FP_mns[,-1])),
                 c(t(Fb_mns[,-1])),c(t(latency_mns[,-1])),c(t(time_mns[,-1])))
all_se1 <- data.frame(scens,rep(methods,4),c(t(TP_sds[,-1])),c(t(FP_sds[,-1])),
                  c(t(Fb_sds[,-1])),c(t(latency_sds[,-1])),c(t(time_sds[,-1])))
colnames(all_mns1) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")
colnames(all_se1) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")


methods2 <- c("CCD","BOCPD","RobOBayes")
scens2 <- c(sapply(5:8,function(i){rep(i,3)}))
all_mns2 <- data.frame(scens2,rep(methods2,4),c(t(TP_mns2[,-1])),c(t(FP_mns2[,-1])),
                 c(t(Fb_mns2[,-1])),c(t(latency_mns2[,-1])),c(t(time_mns2[,-1])))
all_se2 <- data.frame(scens2,rep(methods2,4),c(t(TP_sds2[,-1])),c(t(FP_sds2[,-1])),
                  c(t(Fb_sds2[,-1])),c(t(latency_sds2[,-1])),c(t(time_sds2[,-1])))
colnames(all_mns2) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")
colnames(all_se2) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")

methods3 <- c("CCD","BOCPDMS","rBOCPDMS","BOCPD","RobOBayes")
scens <- c(sapply(9,function(i){rep(i,5)}))
all_mns3 <- data.frame(scens,rep(methods3,1),c(t(TP_mns3[,-1])),c(t(FP_mns3[,-1])),
                       c(t(Fb_mns3[,-1])),c(t(latency_mns3[,-1])),c(t(time_mns3[,-1])))
all_se3 <- data.frame(scens,rep(methods3,4),c(t(TP_sds3[,-1])),c(t(FP_sds3[,-1])),
                      c(t(Fb_sds3[,-1])),c(t(latency_sds3[,-1])),c(t(time_sds3[,-1])))
colnames(all_mns3) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")
colnames(all_se3) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (s)")

all_mns3 <- all_mns3[c(1,4,5),]
all_se3 <- all_se3[c(1,4,5),]


all_mns <- rbind(all_mns1,all_mns2,all_mns3)
all_se <- rbind(all_se1,all_se2,all_se3)

inds<- expand.grid(1:nrow(all_mns),3:ncol(all_mns))

mn_table <- all_mns
sd_table <- all_se
tab <- cbind(all_mns[,1:2],matrix(c(sapply(1:nrow(inds),function(y){
  paste(format(as.numeric(mn_table[inds[y,1],inds[y,2]]),nsmall=2,digits=1),
        "\\textsubscript{",
        format(as.numeric(sd_table[inds[y,1],inds[y,2]]),nsmall=2,digits=1),"}",sep="")
})),nrow(mn_table),ncol(mn_table)-2))

colnames(tab) <- c("Scenario","Method","TP","FP","F-score","Latency","Time (ms)")

tab[,1] <- as.integer(tab[,1])

cap <- "Average metrics (standard errors are in subscripts) for the simulation
study. The scenarios are defined in Table 1.\\\\"
label <- "tab:sim_results"
print(xtable(tab, type = "latex",caption=cap,label=label),
      sanitize.text.function = identity,
      include.rownames = FALSE,
      caption.placement = "top",
      table.placement="H",
      hline.after=c(-1,0,seq(5,20,5),seq(23,35,3)), 
      file = paste("ss_metrics_table_v3.tex",sep=""),
)

make_table(TP_mns,TP_sds,atype="TP",scens=1:4)
make_table(FP_mns,FP_sds,atype="FP",scens=1:4)
make_table(Fb_mns,Fb_sds,atype="F1",scens=1:4)
make_table(latency_mns,latency_sds,atype="latency",scens=1:4)
make_table(time_mns,time_sds,atype="time",scens=1:4)
make_table(TP_mns2,TP_sds2,atype="TP",scens=5:8)
make_table(FP_mns2,FP_sds2,atype="FP",scens=5:8)
make_table(Fb_mns2,Fb_sds2,atype="F1",scens=5:8)
make_table(latency_mns2,latency_sds2,atype="latency",scens=5:8)
make_table(time_mns2,time_sds2,atype="time",scens=5:8)


TP_table <- cbind(1:8,sim_settings,TP_ccdc,TP_bfast,t(TP_bocpd),t(TP_bocpd2),TP_rbocpd)
colnames(TP_table) <- c("Scenario",
                        paste(expression(Delta*mu)),paste(expression(rho)),
                        "b",
                        "CCDC","BOCPDMS","rBOCPDMS","BOCPD","BOCPD-OD")
                     #"BOCPD 1","BOCPD 3","BOCPD 10","BOCPD 50")

FP_table <- cbind(1:8,sim_settings,FP_ccdc,FP_bfast,t(FP_bocpd),t(FP_bocpd2),FP_rbocpd)
colnames(FP_table) <- c("Scenario",
                        paste(expression(Delta*mu)),paste(expression(rho)),
                        "b",
                        "CCDC","BFAST","BOCPD-OD","BOCPD","rBOCPD")
                        #"BOCPD 1","BOCPD 3","BOCPD 10","BOCPD 50")

latency_table <- cbind(1:8,sim_settings,3,latency_bfast,t(latency_bocpd),t(latency_bocpd2))
colnames(latency_table) <- c("Scenario",
                        paste(expression(Delta*mu)),paste(expression(rho)),
                        "b",
                        "CCDC","BFAST","BOCPD-OD","BOCPD")
                       # "BOCPD 1","BOCPD 3","BOCPD 10","BOCPD 50")

print(xtable(TP_table, type = "latex",
             display=c(rep("d",2),rep("f",3),rep("f",4)),digits=c(rep(0,2),1,1,2,rep(3,4))),
      include.rownames = FALSE,
      hline.after=c(-1,0), file = "ss_TP_table.tex")
print(xtable(FP_table, type = "latex",
             display=c(rep("d",2),rep("f",3),rep("f",4)),digits=c(rep(0,2),1,1,2,rep(3,4))),
      include.rownames = FALSE,hline.after=c(-1,0), 
      file = "ss_FP_table.tex")
print(xtable(latency_table, type = "latex",
             display=c(rep("d",2),rep("f",3),rep("f",4)),digits=c(rep(0,2),1,1,2,0,rep(3,3))),
      include.rownames = FALSE,
      hline.after=c(-1,0), file = "ss_latency_table.tex")



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
      if(any(names(cp_info)=="ccdc_cps")){
        if(atype=="ccdc"){
          metrics <- get_metrics_from_cps(cp_info$ccdc_cps,(cp_info$ccdc_cps+2),tol=5)
          t1 <- cp_info$ccdc_time
        }
        if(atype=="bocpdod"){
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

# Plot data generation
npts <- 20
yrs<- 10
t <- seq(1,yrs*365,length.out=npts*yrs)
Y <- generate_scenario(scenario=8,yrs,npts,noerror=F)$Y
fullY <- generate_scenario(scenario=8,yrs,npts,noerror=T)$Y
png("sim_example.png")
par(mfrow=c(2,1),mar=c(2,4,2,2))
plot(Y[,1],ylim=c(0,1),ylab="Response 1",xlab="Time point")
lines(fullY[,1],col="blue")
abline(v=101,col="red")
text(x=101,y=1,labels="101",pos=4)
legend("topright",legend=c("Changepoint","Mean process"),lty=c(1,1),col=c("red","blue"))

plot(Y[,2],ylim=c(0,1),ylab="Response 2",xlab="Time point")
lines(fullY[,2],col="blue")
abline(v=101,col="red")
text(x=101,y=1,labels="101",pos=4)
dev.off()

set.seed(9642)
Y <- c(rnorm(5,0,.5),rnorm(5,2,.5),rnorm(5,-3,.5))
plot(Y,ylim=c(-4,3))
abline(a=-2, b=0)
abline(v=c(6,11))
