load("data/data_base.RDAta")
temp<-data.log.obs%>%pivot_longer(`1`:`28`,names_to="group",values_to = "value")%>%mutate(group=as.numeric(group))
temp<-left_join(temp,sizegroup,by="group")
temps<-sizegroup%>%group_by(group)%>%summarise(count=n())
temp<-left_join(temp,temps,by="group")

ggplot(temp)+geom_boxplot(aes(x=as.numeric(as.character(center)),y=exp(value),group=center,fill=as.factor(count),colour=as.factor(count)),varwidth = TRUE)+
  scale_fill_viridis_d()+scale_colour_viridis_d()+scale_x_log10()+scale_y_log10()+labs(y="Concentration (1/cm^3)",x="Particle size (nm)",fill="Number of \naggregated\nsize bins",colour="Number of \naggregated\nsize bins")

#ggplot(temp)+geom_boxplot(aes(x=size.agg,y=exp(value),group=size.agg))+
#  scale_y_log10()+labs(y="Concentration (1/cm^3)",x="Particle size (nm)")+theme(axis.text.x=element_text(angle=90))

ggsave(file=paste0("plots/data_boxplot.pdf"),device="pdf",width=25,height=15,units="cm")
