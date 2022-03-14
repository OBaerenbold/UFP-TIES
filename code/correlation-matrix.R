library(tidyverse)
data <- read_csv("data/dataset_hr.csv", 
                 col_names = TRUE)

data$index<-1:nrow(data)
head(data)
summary(data)
data.long<-data%>%pivot_longer(cols=`14.6`:`661.2`,names_to = "size",values_to = "value")
#no values below zero
data.long%>%filter(value<0)
#921 values of exactly zero. Fill with half of minimum
data.long%>%filter(value==0)
#join in half of minimum value
data.long.c<-left_join(data.long,data.long%>%filter(value>0)%>%group_by(size)%>%summarise(value_min=min(value,na.rm=T)/2),by="size")
data.long.c<-data.long.c%>%rowwise()%>%mutate(value.c=max(value,value_min))
data.long.c<-data.long.c%>%mutate(size_num=as.numeric(size))
data.long.c.log<-data.long.c%>%mutate(value.c.log=log(value.c))
data.log<-data.long.c.log%>%select(-value,-value_min,-size_num,-value.c)%>%pivot_wider(names_from = size,values_from = value.c.log)
data.lin<-data.long.c.log%>%select(-value,-value_min,-size_num,-value.c.log)%>%pivot_wider(names_from = size,values_from = value.c)

# cor.matrix.log<-as.data.frame(cor(data.log%>%select(-date,-index),use="pairwise.complete.obs"))
# names<-names(data%>%select(-date,-index))
# cor.matrix.log$size1<-names
# cor.matrix.log<-cor.matrix%>%pivot_longer(cols=!c("size1"),names_to="size2",values_to = "cor")

cor.matrix<-as.data.frame(cor(data.lin%>%select(-date,-index),use="pairwise.complete.obs"))
names<-names(data%>%select(-date,-index))
cor.matrix$size1<-names
cor.matrix<-cor.matrix%>%pivot_longer(cols=!c("size1"),names_to="size2",values_to = "cor")

#this part basically does the grouping

low_lim<-0.97 #0.944 for the log one
groups.df<-data.frame(size=names,group=NA)
groups<-list()
z<-1
groups.df$group[1]<-z
groups[[z]]<-names[1]
for(i in 2:length(names)){
  check<-as.logical(cor.matrix%>%filter(size1==groups[[z]][1],size2==names[i])%>%mutate(res= cor>low_lim)%>%pull())
  if(check){
    groups[[z]]<-c(groups[[z]],names[i])
  }
  else{
    z<-z+1
    groups[[z]]<-names[i]
  }
  groups.df$group[i]<-z
}
groups.df<-left_join(groups.df,groups.df%>%group_by(group)%>%summarise(size.agg=ifelse(min(as.numeric(size))==max(as.numeric(size)),size,paste0(min(as.numeric(size))," - ",max(as.numeric(size))))))
center<-round(as.numeric(groups.df%>%group_by(group)%>%summarise(m=mean(as.numeric(size),na.rm=T))%>%pull(m)),1)

cor.matrix<-cor.matrix%>%mutate(gt_lim=as.factor(as.numeric(cor>0.97)))

ggplot(cor.matrix)+geom_tile(aes(x=as.numeric(as.character(size1)),y=as.numeric(as.character(size2)),fill=cor),size=2)+scale_fill_viridis_c()+labs(x="Particle size (nm)",y="Particle size (nm)",fill="Correlation")+scale_x_log10()+scale_y_log10()+
  geom_point(data=groups.df,aes(x=as.numeric(as.character(size)),y=13.5,colour=as.factor(group %% 2)),shape=15,size=2)+scale_colour_discrete(guide="none")

ggplot(cor.matrix)+geom_tile(aes(x=as.numeric(as.character(size1)),y=as.numeric(as.character(size2)),fill=gt_lim))+scale_fill_viridis_d()+labs(x="Particle size (nm)",y="Particle size (nm)",fill="Correlation>0.97")+scale_x_log10()+scale_y_log10()
  


ggsave(file=paste0("plots/correlation-matrix.pdf"),device="pdf",width=20,height=15,units="cm")

ggplot(cor.matrix)+
  geom_tile(aes(x=as.numeric(as.character(size1)),y=as.numeric(as.character(size2)),fill=cor),size=2)+
  labs(x="Particle size (nm)",y="Particle size (nm)",fill="Correlation")+
  scale_x_log10()+scale_y_log10()+ 
  scale_fill_binned(type = "viridis",breaks = c(0,.5,.8,.97,1),guide = guide_coloursteps(even.steps = FALSE),
                    labels=c("<0","[0,0.5)","[0.5,0.8)","[0.8,0.97)","0.97>="))+
  guides(fill = guide_legend(label.position = "right"))+
  geom_point(data=groups.df,aes(x=as.numeric(as.character(size)),y=13.5,colour=as.factor(group %% 2)),shape=15,size=2)+scale_colour_discrete(guide="none")

ggsave(file=paste0("plots/fig1-corr-matrix.pdf"),device="pdf",width=20,height=15,units="cm")
ggsave(file=paste0("plots/fig1-corr-matrix.png"),
       dpi=300,device="png",width=20,height=15,units="cm")
