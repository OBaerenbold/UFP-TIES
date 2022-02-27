library(ggpubr)
library(mcmcplots)
library(nimble)
library(tidyverse)
library(boot)
library(ggmcmc)
#load your data
load("data/data_base.RDAta")
#need matrix only containing concentration, missing data is fine, timestepindex corresponding to row
#called data.log.obs below
#need windspeed and winddirection with no missing - impute with mean at this moment - only small number of missing so shouldn't matter
#data.frame called data.wind with two numerical columns wd and ws, wind direction in degrees
#key for timestepindex, date, hourgroup, etc, called timeref
#key for center of sizegroups called sizegroup
#key for center of hourgroups called hourgroup
#could put all in a single dataframe but no need at this point

#define model
nimble_code <- nimbleCode({
  #order of code is not important
  for(t in 1:Time) {
    for(s in 1:nS){
      # measurement error model
      # mu.Sp is on the logscale here
      # Var are variances 
      species[t,s] ~ dnorm(mu.Sp[t,s], var=Var[s]) 
      
    }
    
    #log_mass[t] is the logarithm of the total concentration at time t
    #inprod is just the matrix multiplication between mu.theta, the source profiles, and
    #p the proportion of the total concentration coming from source s
    for(s in 1:nS){
      mu.Sp[t,s] <-log(inprod(mu.theta[s,1:m],p[1:m,t]))+log_mass[t]
    }
  }
  #each of the source profiles has a non-informative Jeffrey's prior with alpha[i]=1/2. 
  for(j in 1:m){
    mu.theta[1:nS,j]~ddirch(alpha[1:nS])
  }
  
  #A standard 'non-informative' prior on the precision. Easily includes the true values from the data and not of particular interest
  for (s in 1:nS){
    Var[s] ~ dgamma(1, 1)
  }
  
  #some priors
  prec[1]~dgamma(1,1)

  
  #starting point
  for(i in 1:(m-1)){
    #this is the AR(1) process underlying p
    for(t in 1:Time){
    eps[t,i]~dbeta(1,prec[1])
    }
  }
  mean_mass[1]~dnorm(mean=0,sd=10)
  sd_mass[1]~dgamma(1,0.001)
  for(t in 1:(Time)){
    log_mass[t]~dnorm(mean=mean_mass[1],sd=sd_mass[1])
  }
  #Stickbreaking implementation
  for(i in 1:(m-1)){
    knot1[i] ~ dunif(0,max_ws)
    knot2[i] ~ dunif(0, 360)
    invbw1[i] ~ dgamma(1,0.001) # invbw is the inv bandwidth (smoothing parameter)  
    invbw2[i] ~ dgamma(1,0.001) 
  }
  
  for(t in 1:Time){
    #hardcoded transformation from radian to degrees
    for(i in 1:(m-1)){
      #The AR(1) process is transformed to [0,1] and multiplied with the kernel function
      #first kernel function is windspeed, gaussian
      #second kernel function is winddirection, gaussian of sinus^2
      termin[t,i] <- exp(-0.5*invbw1[i]*pow(knot1[i]-ws[t],2))*exp(-0.5*invbw2[i]*pow(sin((knot2[i]-wd[t])*3.14159/360),2))*eps[t,i]
    }
    #stick breaking is a Nimble function that performs the stick-breaking calculation described in the methods section
    p[1:m,t] <- stick_breaking(termin[t,1:(m-1)])
  }
  
})

#the maximum number of sources, start with 10 also try with fewer to test. Doesn't work for 1
maxclust=14
#test with 100 but increase to length of full data, 1604 in our case
maxtime=1604
#which sizegroups should be included in the model. 28 is currently the full after pre-processing
#test with fewer
comps<-1:28

#put all input in a single list that can also be saved if needed
#was done because used to be separate function
model_spec<-list()
#first one is just the code
model_spec[[1]]<-nimble_code
#input data
model_spec[[2]]<-list("species"=as.matrix(data.log.obs)[1:maxtime,comps])


#now all the constants needed in the nimble code above
model_spec[[3]]<-list("Time"=maxtime,"nS"=length(comps), "m"=maxclust,"alpha"=rep(0.5,length(comps)),"wd"=data.wind$wd[1:maxtime], 
               "ws"=data.wind$ws[1:maxtime],"max_ws"=max(data.wind$ws,na.rm=T))
#initial values for some parameters where needed for the mcmc to work. Most are random. This is mostly true for the CAR process
model_spec[[4]]<-list(Var=rep(1,length(comps)),prec=1,log_mass=rep(0,maxtime),eps=matrix(0,nrow=maxtime,ncol=maxclust-1),log_mass=as.vector(matrix(0,nrow=maxtime,ncol=1)))


#the next function first creates the model
#ignore most errors
model <- nimbleModel(model_spec[[1]],
                     data=model_spec[[2]],
                     constants=model_spec[[3]],
                     inits=model_spec[[4]])

## Configure the MCMC 
conf <- configureMCMC(model,useConjugancy=F)
#add monitors to variables of interest
conf$addMonitors("p", "mu.Sp","Var", "mu.theta","log_mass","mean_mass","sd_mass","prec","knot1","knot2","invbw1","invbw2")

## Build the MCMC
mcmc <- buildMCMC(conf)

## Compile the model and the MCMC
cModel <- compileNimble(model)
cMCMC <- compileNimble(mcmc, project = cModel)

## Run the MCMC

samples <- runMCMC(cMCMC, niter=60000, nburnin=30000, thin=30, nchains=1,samplesAsCodaMCMC = T,summary=T)

View(samples$summary)
#single chain at this point to not have to implement lable switching but split the chain to have access to more advanced diagnostics
S_split<-ggs(samples$samples,splitting = T)
S<-ggs(samples$samples)

#this just separates indices from parameter names for easier manipulation
S_prep<-S%>%separate(Parameter,into=c("Parameter","index"),sep="\\[")
S_prep<-S_prep%>%separate(index,into=c("index1","index2"),sep=",")
S_prep<-S_prep%>%separate(index2,into=c("index2","trash"),sep="\\]")%>%select(-trash)
S_prep<-S_prep%>%separate(index1,into=c("index1","trash"),sep="\\]")%>%select(-trash)
S_prep<-S_prep%>%mutate(index1=as.numeric(index1),index2=as.numeric(index2))
#save results before doing too much in case it crashes...
save(samples,maxclust,maxtime,comps,S_prep,file=paste0("samples/","nimble_noar_c",maxclust,"-","_T_",maxtime,"_comps_",min(comps),"-",max(comps),"_",Sys.Date(),".RData"))

#for geweke plots the chain 2 ones are more relevant as they correspond to the second half of the chain in the single chain case
#look at some diagnostics
#mainly traceplots, Rhat, geweke, autocorrelation, crosscorrelation
S_diag<-ggs_diagnostics(S_split)
View(S_diag%>%filter(Diagnostic=="Rhat"))
View(S_diag%>%filter(Diagnostic=="z",Chain==2))
ggs_traceplot(S_split%>%filter(Parameter%in%paste0("p[3, ",c(17,20,9),"]")))
ggs_traceplot(S_split,family="Var")
ggs_geweke(S_split,family="Var")
ggs_geweke(S_split,family="mu.theta")
ggs_traceplot(S_split%>%filter(Parameter=="mu.theta[9, 2]"))
ggs_geweke(S_split,family="knot2")
ggs_geweke(S_split,family="log_mass")
ggs_traceplot(S_split%>%filter(Parameter%in%paste0("log_mass[",5:10,"]")))
ggs_traceplot(S_split%>%filter(Parameter=="sd_mass[1]"))
ggs_traceplot(S_split%>%filter(Parameter%in%paste0("p[1, ",5:10,"]")))
ggs_crosscorrelation(S_split%>%filter(Parameter%in%c(paste0("p[3, ",5:50,"]"))))
ggs_crosscorrelation(S_split%>%filter(Parameter%in%c(paste0("log_mass[",5:50,"]"))))
ggs_crosscorrelation(S_split%>%filter(Parameter%in%c(paste0("p[3, ",1:50,"]"),"psi[3]","mu[1]","prec[1]")))
ggs_traceplot(S_split,family="psi")
