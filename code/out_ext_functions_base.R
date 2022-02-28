clustprob<-function(maxclust,maxtime,S_prep,timeref){
  #aggregate p over iterations to get time-series of proportion
  z<-S_prep%>%filter(Parameter=="p")%>%group_by(index1,index2)%>%summarise(mean=mean(value),low=quantile(value,c(0.025)),high=quantile(value,c(0.975)))
  z<-z%>%rename("source"="index1","time"="index2")
  require(tidyverse)
  #now the concentration needs to be calculated from p and log_mass to create time-series of concentration
  z2<-S_prep%>%filter(Parameter=="log_mass")%>%rename("time"="index1","log_mass"="value")%>%select(-Parameter,-index2)
  z3<-S_prep%>%filter(Parameter=="p")%>%rename("time"="index2","source"="index1","p"="value")%>%select(-Parameter)
  z4<-left_join(z3,z2,by=c("Iteration","time","Chain"))
  z4<-z4%>%mutate(conc=p*exp(log_mass))
  #need the time series of total concentration as output as well
  m<-z4%>%group_by(Iteration,time)%>%summarise(conc=sum(conc))%>%group_by(time)%>%summarise(mean=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  z4<-z4%>%group_by(source,time)%>%summarise(mean=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  #combine concentration and proportion in single dataframe
  z<-left_join(z%>%rename("p_mean"="mean","p_low"="low","p_high"="high"),z4%>%rename("c_mean"="mean","c_low"="low","c_high"="high"),by=c("time","source"))
  
  #create some plots of interest
  plt<-ggplot(z)+geom_line(aes(x=time,y=p_mean,group=source))+geom_ribbon(aes(x=time,ymin=p_low,ymax=p_high,group=source),colour="grey",alpha=0.5)+facet_wrap(~source)
  plt2<-ggplot(z)+geom_line(aes(x=time,y=c_mean,group=source))+geom_ribbon(aes(x=time,ymin=c_low,ymax=c_high,group=source),colour="grey",alpha=0.5)+facet_wrap(~source)+scale_y_log10()
  plt3<-ggplot(m)+geom_line(aes(x=time,y=mean))+geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="grey",alpha=0.5)+scale_y_log10()
  return(list(z,plt,plt2,m,plt3))
}

pollutants_corr<-function(data,poll,maxtime,maxclust,clustfilt){
  library(corrplot)
  data<-data%>%filter(source%in%clustfilt)
  probs<-data%>%pivot_wider(id_cols = c("time","source"),names_from="source",names_prefix = "S_",values_from="p_mean")%>%select(-time)
  cons<-data%>%pivot_wider(id_cols = c("time","source"),names_from="source",names_prefix = "S_",values_from="c_mean")%>%select(-time)
  polls<-poll%>%select(-timeindex)
  polls<-polls[1:maxtime,]
  pdta<-cor(cbind(as.matrix(probs),as.matrix(polls)),use = "pairwise.complete.obs")
  cdta<-cor(cbind(as.matrix(cons),as.matrix(polls)),use = "pairwise.complete.obs")
  levels<-rownames(pdta)
  pdta<-as.data.frame(pdta)
  pdta$type1<-levels
  pdta<-pdta%>%pivot_longer(levels[1]:PMFR,names_to="type2",values_to="cor")
  pdta$type1<-factor(pdta$type1,levels=levels)
  pdta$type2<-factor(pdta$type2,levels=levels)
  clustlength<-length(clustfilt)
  plt1<-ggplot(pdta%>%filter(type1%in%levels[1:clustlength],type2%in%levels[1:clustlength]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred",midpoint=0)+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="black") + 
    ylab('Source')+xlab('Source') + labs(fill="Correlation") + scale_x_discrete(labels=clustfilt) + scale_y_discrete(labels=clustfilt)
  plt2<-ggplot(pdta%>%filter(type1%in%levels[(clustlength+1):length(levels)],type2%in%levels[1:clustlength]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_viridis_c()+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="white") 
  plt3<-ggplot(pdta%>%filter(type1%in%levels[(clustlength+1):length(levels)],type2%in%levels[(clustlength+1):length(levels)]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_viridis_c()+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="white") 
  #same for concentrations
  cdta<-as.data.frame(cdta)
  cdta$type1<-levels
  cdta<-cdta%>%pivot_longer(levels[1]:PMFR,names_to="type2",values_to="cor")
  cdta$type1<-factor(cdta$type1,levels=levels)
  cdta$type2<-factor(cdta$type2,levels=levels)
  plt4<-ggplot(cdta%>%filter(type1%in%levels[1:clustlength],type2%in%levels[1:clustlength]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_viridis_c()+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="white") 
  plt5<-ggplot(cdta%>%filter(type1%in%levels[(clustlength+1):length(levels)],type2%in%levels[1:clustlength]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_viridis_c()+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="white") 
  plt6<-ggplot(cdta%>%filter(type1%in%levels[(clustlength+1):length(levels)],type2%in%levels[(clustlength+1):length(levels)]))+geom_tile(aes(x=type1,y=type2,fill=cor))+scale_fill_viridis_c()+ geom_text(aes(x=type1,y=type2,label = round(cor, 2)),colour="white") 
  
  return(list(list(pdta,plt1,plt2,plt3),list(cdta,plt4,plt5,plt6)))
}

clust_comp<-function(maxclust,maxtime,comps,S_prep,clustfilt=1:maxclust,sizegroup){
  z<-S_prep%>%filter(Parameter=="mu.theta")%>%group_by(index1,index2)%>%summarise(mean=mean(value),low=quantile(value,c(0.025)),high=quantile(value,c(0.975)))
  z<-z%>%rename("source"="index2","size"="index1")
  z<-z%>%mutate(size.center=as.factor(size))
  levels(z$size.center)<-levels(sizegroup$center)[comps]
  z$size.center<-as.numeric(as.character(z$size.center))
  z$source<-as.factor(z$source)
  plt1<-ggplot(z)+geom_pointrange(aes(x=size.center,y=mean,ymin=low,ymax=high))+facet_wrap(~source)
  plt1b<-ggplot(z)+geom_pointrange(aes(x=size.center,y=logit(mean),ymin=logit(low),ymax=logit(high)))+facet_wrap(~source)
  plt1c<-ggplot(z%>%filter(source%in%clustfilt))+geom_line(aes(x=size.center,y=mean,group=source,colour=source))+geom_ribbon(aes(x=size.center,ymin=low,ymax=high,group=source,fill=source),alpha=0.5)+
    scale_colour_viridis_d()+scale_x_log10()+xlab('Particle Size (nm)')+ ylab('Proportion of particle conc.')+ggtitle('Source profile') +labs(fill='Source',colour='Source')
  z2<-S_prep%>%filter(Parameter%in%c("invbw1","invbw2","knot1","knot2"))%>%group_by(index1,Parameter)%>%summarise(mean=mean(value),low=quantile(value,c(0.025)),high=quantile(value,c(0.975)))
  z2<-z2%>%rename("source"="index1")
  z3<-S_prep%>%filter(Parameter%in%c("invbw1","invbw2","knot1","knot2"))%>%select(-index2,-Chain)
  z3<-z3%>%rename("source"="index1")
  z3<-z3%>%pivot_wider(id_cols=c("Iteration","Parameter","source"),names_from=Parameter,values_from=value)
  #only pick 100 samples to predict mean surface
  samp<-round(seq(from=100,to=(max(S_prep$Iteration)-100),len=100))
  z3<-z3%>%filter(Iteration%in%samp)
  #create grid for wind coordinates prediction
  wsgrid<-seq(0,6,by=0.4)
  wdgrid<-seq(0,359.9,by=10)
  wswd<-data.frame(expand.grid(wsgrid,wdgrid))
  names(wswd)<-c("wsgrid","wdgrid")
  #merge grid into z3
  z3<-merge(z3,wswd,by=NULL)
  z3<-z3%>%mutate(kernel=exp(-0.5*invbw1*(knot1-wsgrid)^2)*exp(-0.5*invbw2*(sin((knot2-wdgrid)*3.14159/360))^2))
  z3<-z3%>%group_by(wsgrid,wdgrid,source)%>%summarise(kernel=mean(kernel))
  plt2<-ggplot(z3%>%filter(source%in%clustfilt))+geom_tile(aes(x=wdgrid,y=wsgrid,fill=kernel))+scale_fill_viridis_c(option = "magma",limits=c(0,1))+
    coord_polar(start = -5 / 180 * pi) +geom_point(aes(x=0,y=0),size=2,colour="red")+facet_wrap(~source)+ylab('Wind Speed (m/s)')+xlab('Wind direction')+labs(fill='Kernel')+
    scale_x_continuous(breaks = c(0,90,180,270),labels = c("N", "E", "S","W"))  
  
  return(list(list(z,plt1,plt1b,plt1c),list(z3,plt2)))
  }

AggregateMeans<-function(maxclust,maxtime,S_prep,timeref){
  z1<-S_prep%>%filter(Parameter=="log_mass")%>%rename("time"="index1","log_mass"="value")%>%select(-Parameter,-index2)
  z2<-S_prep%>%filter(Parameter=="p")%>%rename("time"="index2","source"="index1","p"="value")%>%select(-Parameter)
  z3<-left_join(z2,z1,by=c("Iteration","time","Chain"))
  z4<-z3%>%mutate(conc=p*exp(log_mass))
  z4<-left_join(z4,timeref%>%mutate(time=timeindex),by="time")
  #calculate summaries by cluster
  a1<-z4%>%group_by(source,Iteration,weekday)%>%summarise(conc=mean(conc))%>%
    group_by(weekday,source)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  a2<-z4%>%group_by(source,Iteration,month)%>%summarise(conc=mean(conc))%>%
    group_by(month,source)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  a3<-z4%>%group_by(source,Iteration,center_hour)%>%summarise(conc=mean(conc))%>%
    group_by(center_hour,source)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  #just calucate the totals as well
  b1<-z4%>%group_by(Iteration,timeindex,weekday)%>%summarise(conc=sum(conc))%>%group_by(Iteration,weekday)%>%summarise(conc=mean(conc))%>%
    group_by(weekday)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  b2<-z4%>%group_by(Iteration,timeindex,month)%>%summarise(conc=sum(conc))%>%group_by(Iteration,month)%>%summarise(conc=mean(conc))%>%
    group_by(month)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  b3<-z4%>%group_by(Iteration,timeindex,center_hour)%>%summarise(conc=sum(conc))%>%group_by(Iteration,center_hour)%>%summarise(conc=mean(conc))%>%
    group_by(center_hour)%>%summarise(m=mean(conc),low=quantile(conc,c(0.025)),high=quantile(conc,c(0.975)))
  b1$source=0
  b2$source=0
  b3$source=0
  c1<-rbind(a1,b1)
  c2<-rbind(a2,b2)
  c3<-rbind(a3,b3)
  c1$source_n<-as.factor(c1$source)
  levels(c1$source_n)<-c("Total",paste0("Source ",1:maxclust))
  c2$source_n<-as.factor(c2$source)
  levels(c2$source_n)<-c("Total",paste0("Source ",1:maxclust))
  c3$source_n<-as.factor(c3$source)
  levels(c3$source_n)<-c("Total",paste0("Source ",1:maxclust))
  plt1<-ggplot(c1)+geom_pointrange(aes(x=weekday,y=m,ymin=low,ymax=high,group=source_n))+scale_y_log10()+facet_wrap(~source_n,scales="free_y")
  plt2<-ggplot(c2)+geom_pointrange(aes(x=month,y=m,ymin=low,ymax=high,group=source_n))+scale_y_log10()+facet_wrap(~source_n,scales="free_y")
  
  plt3<-ggplot(c3)+geom_pointrange(aes(x=center_hour,y=m,ymin=low,ymax=high,group=source_n))+scale_y_log10()+facet_wrap(~source_n,scales="free_y")
  return(list(c1,c2,c3,plt1,plt2,plt3))
}

library(scales)
logit_perc <- trans_new("logit perc",
                        transform = function(x)qlogis(x/100),
                        inverse = function(x)100*plogis(x)
)

Comp_summary<-function(maxclust,agg_means,cl_comp,clustprobs,pollcorr){
  z<-cl_comp[[1]][[1]]
  z3<-cl_comp[[2]][[1]]
  c1<-agg_means[[1]]
  c2<-agg_means[[2]]
  c3<-agg_means[[3]]
    f<-clustprobs[[1]]
    fs1<-f%>%group_by(source)%>%summarise(c=sum(c_mean))
    fs1<-fs1%>%mutate(c_p=c/sum(fs1$c))
  pcor<-pollcorr[[1]][[1]]
  ccor<-pollcorr[[2]][[1]]
  t<-list()
  for(i in 1:maxclust){
  cor<-rbind(pcor%>%filter(type1==paste0("S_",i))%>%select(type2,cor)%>%mutate(type="prop."),ccor%>%filter(type1==paste0("S_",i))%>%select(type2,cor)%>%mutate(type="conc."))
  cor<-cor%>%separate(type2,into=c("a","b"),sep=1)%>%filter(a!="S")%>%mutate(type2=paste0(a,b))%>%select(-a,-b)
  pltcor<-ggplot(cor)+geom_point(aes(x=type2,y=cor,group=type,colour=type,shape=type))+labs(x="",y="Correlation",shape="",colour="")+ theme(axis.text.x = element_text(angle = 90))
  max<-round(max(c(c1%>%filter(source==i)%>%pull(high),c2%>%filter(source==i)%>%pull(high),c3%>%filter(source==i)%>%pull(high)))*1.1,-1)
  min<-round(min(c(c1%>%filter(source==i)%>%pull(low),c2%>%filter(source==i)%>%pull(low),c3%>%filter(source==i)%>%pull(low)))*0.9,-1)
  #prop<-ggplot(f%>%filter(clust==i))+geom_line(aes(x=time,y=p_mean))+geom_ribbon(aes(x=time,ymin=p_low,ymax=p_high),colour="grey",alpha=0.5)
  #conc<-ggplot(f%>%filter(clust==i))+geom_line(aes(x=time,y=c_mean))+geom_ribbon(aes(x=time,ymin=c_low,ymax=c_high),colour="grey",alpha=0.5)+scale_y_log10()
  prof<-ggplot(z%>%filter(source==i))+geom_pointrange(aes(x=size.center,y=mean,ymin=low,ymax=high))+scale_x_log10()+labs(x="Particle Size (nm)",y='Profile')
  prof_log<-ggplot(z%>%filter(source==i))+geom_pointrange(aes(x=size.center,y=mean,ymin=low,ymax=high))+ coord_trans(y = logit_perc)+scale_y_continuous(breaks=c(0.0001,0.001,0.01,0.5,0.1,0.2))+scale_x_log10()+labs(x="Particle Size (nm)",y='Logit profile')
  if(i<maxclust){
  wind<-ggplot(z3%>%filter(source==i))+geom_tile(aes(x=wdgrid,y=wsgrid,fill=kernel))+scale_fill_viridis_c(option = "magma",limits=c(0,1))+
    coord_polar(start = -5 / 180 * pi) +geom_point(aes(x=0,y=0),size=2,colour="red")+ylab('Wind Speed (m/s)')+xlab('Wind direction')+
    labs(fill='Kernel')+scale_x_continuous(breaks = c(0,90,180,270),labels = c("N", "E", "S","W")) 
  }
  week<-ggplot(c1%>%filter(source==i))+geom_pointrange(aes(x=weekday,y=m,ymin=low,ymax=high))+scale_y_log10(limits=c(min,max))+labs(x="Weekday",y="Mean conc. (1/cm^3)")
  month<-ggplot(c2%>%filter(source==i))+geom_pointrange(aes(x=month,y=m,ymin=low,ymax=high))+scale_y_log10(limits=c(min,max))+labs(x="Month",y="Mean conc. (1/cm^3)")+ theme(axis.text.x = element_text(angle = 90))
  day<-ggplot(c3%>%filter(source==i))+geom_pointrange(aes(x=center_hour,y=m,ymin=low,ymax=high))+scale_y_log10(limits=c(min,max))+xlim(0,23)+labs(x="Time of day (h)",y="Mean conc. (1/cm^3)")
  if(i<maxclust){
  t[[i]]<-annotate_figure(ggarrange(   # plot4 in first row
    ggarrange(wind,day, week,month, ncol = 4),
    ggarrange(prof,prof_log,pltcor,ncol=3),
    nrow = 2  # plot1 and plot2 in second row
  ),top=paste0("Source-",i," / ",round(fs1$c_p[i]*100,1),"% of total concentration"))
  }else{
    t[[maxclust]]<-annotate_figure(ggarrange(   # plot4 in first row
      ggarrange(day, week,month, ncol = 3),
      ggarrange(prof,prof_log,pltcor,ncol=3),
      nrow = 2  # plot1 and plot2 in second row
    ),top=paste0("Source-",i," / ",round(fs1$c_p[i]*100,1),"% of total concentration"))
  }
  }
  return(t)
}

Total_summary<-function(maxclust,agg_means,clustprobs,S_prep,sizegroup,timeref){
  v<-S_prep%>%filter(Parameter=="Var")%>%group_by(index1)%>%summarise(mean=mean(value),low=quantile(value,c(0.025)),high=quantile(value,c(0.975)))
  v<-left_join(v,unique(sizegroup%>%select(index1=group,center)),by="index1")
  v<-v%>%rename("size"="index1")
  v$center<-as.numeric(as.character(v$center))

  c1<-agg_means[[1]]
  c2<-agg_means[[2]]
  c3<-agg_means[[3]]
  f1<-clustprobs[[1]]
  f1<-left_join(f1%>%ungroup(),timeref%>%ungroup()%>%select(time=timeindex,center_hour,weekday,month),by="time")
  ft1<-f1%>%group_by(time,center_hour,weekday,month)%>%summarise(c_mean=sum(c_mean))
  fs1<-f1%>%group_by(source)%>%summarise(c=sum(c_mean))
  fs1<-fs1%>%mutate(c_p=c/sum(fs1$c))
  #f2<-clustprobs[[4]]
  week<-ggplot(ft1)+geom_boxplot(aes(x=weekday,y=c_mean))+scale_y_log10(limits=c(10,70000))+geom_pointrange(data=c1%>%filter(source_n=="Total"),aes(x=weekday,y=m,ymin=low,ymax=high),colour="red",size=0.3)+labs(x="Weekday",y="c(t)")+ theme(axis.title.y = element_text(angle=0))
  month<-ggplot(ft1)+geom_boxplot(aes(x=month,y=c_mean))+scale_y_log10(limits=c(10,70000))+geom_pointrange(data=c2%>%filter(source_n=="Total"),aes(x=month,y=m,ymin=low,ymax=high),colour="red",size=0.3)+labs(x="Month",y="c(t)")+ theme(axis.title.y = element_text(angle=0))
  day<-ggplot(ft1)+geom_boxplot(aes(x=center_hour,y=c_mean,group=center_hour))+scale_y_log10(limits=c(10,70000))+geom_pointrange(data=c3%>%filter(source_n=="Total"),aes(x=center_hour,y=m,ymin=low,ymax=high),colour="red",size=0.3)+xlim(0,23)+labs(x="Time of day (h)",y="c(t)")+ theme(axis.title.y = element_text(angle=0))
  #barp<-ggplot(f1, aes(x=as.factor(clust),y=p)) + 
  #  geom_bar(stat="identity")+labs(y="Mean prop. of conc.",x="Cluster")
  prop<-ggplot(fs1)+geom_point(aes(x=factor(source),y=c_p*100))+geom_hline(aes(yintercept=1),colour="red")+
    labs(x="Source",y="Prop. of conc. (%)")+scale_y_continuous(breaks=c(0,1,5,10,15,20,25,30,40,50))
  barc<-ggplot(f1, aes(x=as.factor(source),y=c)) + 
    geom_boxplot(aes(group=source,y=c_mean))+labs(y="c(t)*f(k,t)",x="Source") 
  #conc<-ggplot(f2)+geom_line(aes(x=time,y=mean))+geom_ribbon(aes(x=time,ymin=low,ymax=high),colour="grey",alpha=0.5)+scale_y_log10()
  sds<-ggplot(v)+geom_pointrange(aes(x=center,y=sqrt(mean),ymin=sqrt(low),ymax=sqrt(high)))+labs(y=expression(sigma[p]),x="Size bin center (nm)")+scale_x_log10()
  plt<-ggarrange(   # plot4 in first row
    ggarrange(prop,barc,sds,ncol=3),
    annotate_figure(ggarrange(day,week,month, ncol = 3),
                    top = text_grob("Aggregate mean", size = 12,hjust = 2.5)),
    nrow = 2  # plot1 and plot2 in second row
  )
  return(plt)
}

check_fit_bysize<-function(maxclust,maxtime,comps,S_prep,data.log.obs){
  z<-S_prep%>%filter(Parameter=="mu.Sp")%>%group_by(index1,index2)%>%summarise(mean=mean(value),low=quantile(value,c(0.025)),high=quantile(value,c(0.975)))
  z<-z%>%rename("time"="index1","size"="index2")
  t<-list()
  n<-1
  for(i in comps){
    data<-data.frame(value=data.log.obs[1:maxtime,i],time=1:maxtime)
    names(data)[1]="value"
    data<-left_join(data,z%>%filter(size==i))
    t[[n]]<-ggplot(data)+geom_point(aes(x=time,y=value))+geom_line(aes(x=time,y=mean),colour="red")+geom_ribbon(aes(x=time,ymin=low,ymax=high),fill="red",alpha=0.5)+labs(title=paste0("Size ",levels(sizegroup$size.agg)[i]," - on log-scale"))
    n<-n+1
  }
  return(t)
}
