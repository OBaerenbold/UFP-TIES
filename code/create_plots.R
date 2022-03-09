library(ggpubr)
library(mcmcplots)
library(nimble)
library(tidyverse)
library(boot)
library(ggmcmc)
#load your data
load("data/data_base.RDAta")
#load the samples
#name<-"nimble_noar_a2_c14-_T_1604_comps_1-28_2022-03-01"
#name<-"nimble_nowind_r1_c14-_T_1604_comps_1-28_2022-03-07"
name<-"nimble_noar_r1_c14-_T_1604_comps_1-28_2022-03-08"
load(paste0("samples/",name,".RData"))
dir.create(paste0("plots/",name))
#load some functions to extract results
source(file = "code/out_ext_functions_base.R")
#clustfilt<-c(1,3,4,7,9,10,12,13,14) a2
#clustfilt<-c(1:4,6,9,10,13,14)
clustfilt<-c(1,7,9,10,11,12,13,14)
#names<-c("Airport","Urban","Sec. aerosols B","UKN winter","Fresh traffic","UKN 1","Sec. aerosols A","Aged traffic","UKN 2")
names<-as.character(clustfilt)
#some of the output of functions is used as input for the next ones so care is needed when making chances to output of functions

#creating time-series of proportion of concentration and concentration of each cluster
clustprobs=clustprob(S_prep=S_prep,maxtime=maxtime,maxclust = maxclust)
clustprobs[[2]]
ggsave(file=paste0("plots/",name,"/prop-timeseries",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/prop-timeseries",".pdf"),device="pdf",width=20,height=15,units="cm")
clustprobs[[3]]
ggsave(file=paste0("plots/",name,"/conc-timeseries",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/conc-timeseries",".pdf"),device="pdf",width=20,height=15,units="cm")
clustprobs[[5]]
ggsave(file=paste0("plots/",name,"/total-conc-timeseries",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/total-conc-timeseries",".pdf"),device="pdf",width=20,height=15,units="cm")
#calculates correlation matrix of time series of each cluster with pollutations
pollcorr<-pollutants_corr(data=clustprobs[[1]],poll=poll,maxtime=maxtime,maxclust=maxclust,clustfilt=clustfilt)
pollcorr[[1]][[2]]
ggsave(file=paste0("plots/",name,"/clust-prob-corr",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clust-prob-corr",".pdf"),device="pdf",width=20,height=15,units="cm")
pollcorr[[1]][[3]]
ggsave(file=paste0("plots/",name,"/clust-prob-poll-corr",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clust-prob-poll-corr",".pdf"),device="pdf",width=20,height=15,units="cm")
pollcorr[[2]][[2]]
ggsave(file=paste0("plots/",name,"/clust-conc-corr",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clust-conc-corr",".pdf"),device="pdf",width=20,height=15,units="cm")
pollcorr[[2]][[3]]
ggsave(file=paste0("plots/",name,"/clust-conc-poll-corr",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clust-conc-poll-corr",".pdf"),device="pdf",width=20,height=15,units="cm")

pollcorr[[1]][[4]]
ggsave(file=paste0("plots/",name,"/poll-corr",".png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/poll-corr",".pdf"),device="pdf",width=20,height=15,units="cm")

#calculates combined plots of cluster profiles, need to select which ones in plot with clustfilt

#cl_comp<-clust_comp(maxclust,maxtime,comps,S_prep,clustfilt=clustfilt,sizegroup,names=names)
#clust_comp_gen as in generic which is not adapted to the specific run we use in the results
cl_comp<-clust_comp_gen(maxclust,maxtime,comps,S_prep,clustfilt=clustfilt,sizegroup,names=names)
#clust_comp_nowind for the runs without wind data
#cl_comp<-clust_comp_nowind(maxclust,maxtime,comps,S_prep,clustfilt=clustfilt,sizegroup,names=names)

cl_comp[[1]][[4]]
ggsave(file=paste0("plots/",name,"/clusts-profile.png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clusts-profile.pdf"),device="pdf",width=20,height=15,units="cm")

cl_comp[[2]][[2]]
ggsave(file=paste0("plots/",name,"/clusts-windkernel.png"),device="png",width=20,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/clusts-windkernel.pdf"),device="pdf",width=20,height=15,units="cm")

#calculates combined plots of aggregated concentration over weekday, month, hour
agg_means<-AggregateMeans(maxclust,maxtime,S_prep,timeref)
agg_means[[4]]
ggsave(file=paste0("plots/",name,"/agg-weekday.png"),device="png",width=25,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/agg-weekday.pdf"),device="pdf",width=25,height=15,units="cm")
agg_means[[5]]
ggsave(file=paste0("plots/",name,"/agg-month.png"),device="png",width=25,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/agg-month.pdf"),device="pdf",width=25,height=15,units="cm")
agg_means[[6]]
ggsave(file=paste0("plots/",name,"/agg-day.png"),device="png",width=25,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/agg-day.pdf"),device="pdf",width=25,height=15,units="cm")

#creates summary plots by cluster
comps_sum<-Comp_summary(maxclust,agg_means,cl_comp,clustprobs,pollcorr)
#just in case of results without wind kernel
#comps_sum<-Comp_summary_nowind(maxclust,agg_means,cl_comp,clustprobs,pollcorr)
for(i in 1:maxclust){
  print(comps_sum[[i]])
  ggsave(file=paste0("plots/",name,"/clust-profile-",i,".png"),device="png",width=25,height=15,units="cm")
  ggsave(file=paste0("plots/",name,"/clust-profile-",i,".pdf"),device="pdf",width=25,height=15,units="cm")
}

#creates overview plot
tot_sum<-Total_summary(maxclust,agg_means,clustprobs,S_prep,sizegroup,timeref,hourgroup)
tot_sum
ggsave(file=paste0("plots/",name,"/total-summary.png"),device="png",width=25,height=15,units="cm")
ggsave(file=paste0("plots/",name,"/total-summary.pdf"),device="pdf",width=25,height=15,units="cm")

#Compares prediction for each sizegroup with data
fit_sizes<-check_fit_bysize(maxclust,maxtime,comps,S_prep,data.log.obs)
for(i in comps){
  print(fit_sizes[[i]])
  ggsave(file=paste0("plots/",name,"/sizegroup-",i,".png"),device="png",width=25,height=15,units="cm")
  ggsave(file=paste0("plots/",name,"/sizegroup-",i,".pdf"),device="pdf",width=25,height=15,units="cm")
}

