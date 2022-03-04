library(data.table)
library(roboBayes)
library(dplyr)

####### Figure 1: Time series of SWIR2 and NDVI with filtered points shown
# visualization data (point 70)
# Myanmar data
### with bad points
# Load visualization dataset
# variables needed ####
set.seed(42)
country <- "myanmar" #all scripts
dataFolder <- paste0("data/", country, "/") #all scripts
inputFolder <- "visualization_points" 
trainingPointsFile <- "visualization_points.csv"

#load pointIDs
base <- fread(paste0(dataFolder, trainingPointsFile))
pointids <- sort(unique(base[,pointid]))

pointData <- trainingPointsFile

fn <- paste0(dataFolder,inputFolder,"/landsatsr.csv")
df <- read.table(fn,sep=",",header=T)
df$date <- as.Date(df$date,format="%Y-%m-%d")

#update dates
date_fn <- function(dae){
  date_true <- ifelse(is.na(as.Date(dae)) | (Sys.Date() - as.Date(dae)) > 10000, 
                      as.Date(dae, format="%d/%m/%Y"),
                      as.Date(dae))
  date_true <- as.Date(date_true, origin="1970-01-01")
  return(date_true)
}

base <- base[, `:=` (date_dist = date_fn(dateDist))]

get_moddate <- function(sat_list=sat_list){
  dat <- sat_list
  setkey(dat, pointid)
  setkey(base, pointid)
  dat <- dat[base]
}

df <- get_moddate(data.table(df))

df$B11 <- df$B11*0.0000275-0.2

ndvi <- transmute(df,ndvi=(B5-B4)/(B5+B4))
swir2 <- transmute(df,swir2=B7*1e-4)
goodquality <- transmute(df,goodquality=((pixel_qa==322)&(radsat_qa==0)&(sr_aerosol<=96)))

df <- cbind(df,ndvi,swir2,goodquality)
dfr <- dplyr::select(df,c("pointid","date_dist","date","goodquality","ndvi","swir2"))
#dfr <- dplyr::filter(dfr,valid_cloudAndRadqa==T)
# split into list of dfs
dfc <- split(dfr,by="pointid")

# split into point ids for pure dataset
#df_pure <- df %>% mutate(goodquality=valid)
#4. filter out bad qa
df_pure <- df_pure %>% filter(goodquality)

dfl_pure <- split(df_pure,df_pure$pointid)
dfl_pure <- lapply(dfl_pure,function(x){
  xi <- mutate(x,t=as.numeric(x$date-min(x$date)))
  xi <- arrange(xi,date)
  mutate(xi,ti=1:nrow(xi))
  # filter
})
#3. Keep in some bad points
corrupted <- sample(which(!df$goodquality&!is.na(df$ndvi)&!is.na(df$swir2)),200)
#df <- df %>% mutate(goodquality=valid)
df$goodquality[corrupted] <- T
#4. filter out bad qa
df <- df %>% filter(goodquality)

#df <- df %>% filter(valid_qa&valid_radsatqa&valid_aerosol&valid_cloudAndRadqa)

#4. split into pointids and sort
dfl <- split(df,df$pointid)
dfl <- lapply(dfl,function(x){
  xi <- mutate(x,t=as.numeric(x$date-min(x$date)))
  xi <- arrange(xi,date)
  mutate(xi,ti=1:nrow(xi))
  # filter
})

# ### Figure 1: SWIR, NDVI ts with filtered points
# #2. calculate NDVI and SWIR
# dfc <- dfc
# dfc <- F1_onboard$sat_list$l8sr
# dfc <- dfc %>% mutate(ndvi = (B5-B4)/(B5+B4))
# dfc <- dfc %>% mutate(swir2 = B7/1e4)
# dfc <- dfc %>% mutate(goodquality=valid_qa&valid_radsatqa&valid_aerosol&valid_cloudAndRadqa)
# 
# dfcl <- split(dfc,dfc$pointid)
#95 is okay, 92, 70, 68,65, 54. 70 is  great!
dfcl <- dfc
pit <- 70
dfi <- dfcl[[pit]]
png("plots/Myanmar_data_swir2_QA.png",width=500,height=200)
par(mfrow=c(1,1),mar=c(2,4,1,1),cex=1.5)
plot(dfi$date,dfi$swir2,ylim=c(0,1),xlab="Year",ylab="SWIR2")
points(dfi$date[dfi$goodquality==F],dfi$swir2[dfi$goodquality==F],pch=4,col=2,lwd=2)
points(dfi$date[44],dfi$swir2[44],pch=1,cex=2,col=4,lwd=2)

dev.off()

png("plots/Myanmar_data_ndvi_QA.png",width=500,height=200)
par(mfrow=c(1,1),mar=c(2,4,1,1),cex=1.5)
plot(dfi$date,dfi$ndvi,ylim=c(0,1),xlab="Year",ylab="NDVI")
points(dfi$date[dfi$goodquality==F],dfi$ndvi[dfi$goodquality==F],pch=4,col=2,lwd=2)
points(dfi$date[44],dfi$ndvi[44],pch=1,cex=2,col=4,lwd=2)
dev.off()

png("plots/Myanmar_legend_QA.png",width=150,height=100)
par(mfrow=c(1,1),mar=c(0,0,0,0),cex=1.5)
plot.new()
legend("topleft",xpd=T,legend=c("Data","Invalid QA","Outlier"),col=c(1,2,4),pch=c(1,4,1),lwd=c(1,2,2),lty=c(NA,NA,NA))
dev.off()


####### Figure 2: Observed SWIR2 and NDVI measurements
# visualization data (point 24)
pit <- 24
dfi <- dfl[[pit]]

png("Plots/Myanmar_data.png",width=1000,height=200)
par(mfrow=c(1,2),mar=c(2,4,1,1),cex=1.5)
plot(dfi$date,dfi$swir2,ylim=c(0,1),xlab="Year",ylab="SWIR2")
abline(v=dfi$date_dist,col=2)
legend("topleft",legend=c("Data","Disturbance"),bty="n",lty=c(NA,1),col=c(1,2),pch=c(1,NA))
plot(dfi$date,dfi$ndvi,ylim=c(0,1),xlab="Year",ylab="NDVI")
abline(v=dfi$date_dist,col=2)
dev.off()



####### Figure 3: Simulation signals with an outlier, run length plot, 
# and state diagram

set.seed(9264)
source("code/functions/sim_study_helpers.R")
scenario <- 7
# generate data
yrs <- 15
npts=18
t <- seq(1,yrs*365,length.out=yrs*npts)
period=365

X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),t/1000)

make_scenario <-  generate_scenario(scenario,yrs=yrs,npts=npts,l1=10,l2=5)

Y <- make_scenario$Y

rl <- c(seq(1,180),seq(1,90))
st <- c(rep(1,180),rep(2,90))
st[make_scenario$outlier] <- 0

png("plots/simulation_data_7_signal1.png",width=4,height=1.7,units="in",res=72)
par(mfrow=c(1,1),mar=c(4,4,0.5,1.7),cex=1)
# make labels and margins smaller
plot(1:270,Y[,1],ylab="Signal 1",xaxt="n",xlab="",pch=20)
text(make_scenario$outlier,0.78,"outlier",pos=4)
text(300,0.8,"(a)",xpd=NA)
dev.off()

png("plots/simulation_data_7_state.png",width=4,height=1.7,units="in",res=72)
par(mfrow=c(1,1),mar=c(4,4,0.5,1.7),cex=1)
plot(1:270,st,type="l",ylim=c(0,3),yaxt='n',ylab="State",xaxt="n",xlab="")
axis(2,at=c(0,1,2))
text(300,3,"(b)",xpd=NA)
dev.off()

png("plots/simulation_data_7_signal2.png",width=4,height=1.7,units="in",res=72)
par(mfrow=c(1,1),mar=c(4,4,0.5,1.7),cex=1)
plot(1:270,Y[,2],ylab="Signal 2",xlab="Time Point",pch=20)
text(make_scenario$outlier,0.14,"outlier",pos=4)
text(300,0.65,"(c)",xpd=NA)
dev.off()

png("plots/simulation_data_7_rl.png",width=4,height=1.7,units="in",res=72)
par(mfrow=c(1,1),mar=c(4,4,0.5,1.7),cex=1)
plot(rl,type='l',ylab="Run length",xlab="Time Point")
text(300,180,"(d)",xpd=NA)

dev.off()


####### Figure 4: workflow

####### Figure 5: Myanmar results: (visualization data point 70)
# with outliers removed and change dates shown
load("code/paper plots/myanmar_default.RData")
load("code/paper plots/myanmar_priors.RData")
pit <- 1#73
dfi <- dfl[[pit]]
t <- dfi$t
period <- 365
# run bocpdod and pull out R
X <- cbind(1,sin(2*pi*t/period),cos(2*pi*t/period),(t-min(t))/(max(t)-min(t)))

Yi <- cbind(dfi$ndvi,dfi$swir2)

# pull out change point number
truecp <- dfi$ti[which(dfi$date_dist<=dfi$date)[1]]
piMine <- priors


bocpd_mod_o <- roboBayes::roboBayes(datapts = Yi,
                     covariates = X,
                     RoboBayes = NULL,
                     par_inits = piMine,
                     Lsearch = 15,
                     Lwindow = 8,
                     Lgroup=5,
                     lambda=100,
                     cpthresh = 0.8,
                     truncRmin = 300,
                     cptimemin = 12,
                     Lm = 8,
                     cp_delay = 3,
                     kt = 1,
                     pc=0.5,
                     alpha=0.9,
                     getR = TRUE,
                     getOutliers = TRUE,
                     getModels = FALSE)



# analyze using BOCPD
bocpd_mod_no <- roboBayes::roboBayes(datapts = Yi,
                                    covariates = X,
                                    RoboBayes = NULL,
                                    par_inits = piMine,
                                    Lsearch = 15,
                                    Lwindow = 8,
                                    Lgroup=5,
                                    lambda=100,
                                    cpthresh = 0.8,
                                    truncRmin = 300,
                                    cptimemin = 12,
                                    Lm = 8,
                                    cp_delay = 3,
                                    kt = 1,
                                    pc=0.5,
                                    alpha=0.9,
                                    getR = TRUE,
                                    getOutliers = FALSE,
                                    getModels = FALSE)


# put outlier results in dfir
dfio <- dfi[bocpd_mod_o$outliers,]

library(fields)
bocpd_mod <- bocpd_mod_o
s1=expand.grid(1:nrow(bocpd_mod$RFull),1:ncol(bocpd_mod$RFull))
thing1 <- data.frame(s1,R2=c(t(bocpd_mod$RFull)))
names(thing1) <- c("s1","s2","R2")
bocpd_mod <- bocpd_mod_no
s1=expand.grid(1:nrow(bocpd_mod$RFull),1:ncol(bocpd_mod$RFull))
thing2 <- data.frame(s1,R2=c(t(bocpd_mod$RFull)))
names(thing2) <- c("s1","s2","R2")

dfi$outliers_ndvi <- NA
dfi$outliers_ndvi[bocpd_mod$outliers] <- dfi$ndvi[bocpd_mod$outliers]

library(ggplot2)

# now plot clean data for comparison
bocpd_mod <- bocpd_mod_o
png("Plots/Myanmar_data_1_swir2.png",width=4,height=1.5,units="in",res=1200)
par(mfrow=c(1,1),mar=c(4,4,0.1,0.1),cex=0.7)
plot(dfi$date,dfi$swir2,ylim=c(0,1),xlab="Year",ylab="SWIR2")
abline(v=dfi$date[unique(c(bocpd_mod$cpInds))[-1]],col="red",lwd=2)
points(dfi$date[bocpd_mod$outliers],dfi$swir2[bocpd_mod$outliers],pch=4,col="red",lwd=2)
legend("topleft",legend=c("data","outliers","disturbance"),bty="n",lty=c(NA,NA,1),lwd=c(NA,2,2),col=c(1,"red","red"),pch=c(1,4,NA))
dev.off()
png("Plots/Myanmar_data_1_ndvi.png",width=4,height=1.5,units="in",res=1200)
par(mfrow=c(1,1),mar=c(4,4,0.1,0.1),cex=0.7)
plot(dfi$date,dfi$ndvi,ylim=c(0,1),xlab="Year",ylab="NDVI")
abline(v=dfi$date[unique(c(bocpd_mod$cpInds))[-1]],col="red",lwd=2)
points(dfi$date[bocpd_mod$outliers],dfi$ndvi[bocpd_mod$outliers],pch=4,col="red",lwd=2)
dev.off()



####### Figure 6: Myanmar run lengths
######## Do this one twice, on with outliers and one without
plt <- ggplot(thing2)
plt <- plt + theme_classic() 
plt <- plt + scale_color_gradient(low="white", high="black",
                                  name = "Run length \n probability \n")
plt <- plt + geom_point(aes(x=s1,y=s2,color=R2))
plt <- plt + xlab("Time point") + ylab("Run length")
plt
ggsave(filename="plots/myanmar_rl.png",width=3,height=1,units="in",scale=1.2)

plt <- ggplot(thing1)
plt <- plt + theme_classic() 
plt <- plt + scale_color_gradient(low="white", high="black",
                                  name = "Run length \n probability \n")
plt <- plt + geom_point(aes(x=s1,y=s2,color=R2))
plt <- plt + xlab("Time point") + ylab("Run length")
plt
ggsave(filename="plots/myanmar_rl_outliers.png",width=3,height=1,units="in",scale=1.2)
