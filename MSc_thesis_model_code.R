#Mathematical model for competition growth
rm(list=ls())
setwd("C:/model")
# 0.Load packages ####
library(pacman)
p_load(deSolve, tidyverse)
library(readxl)


# 1.Load data ####
## absenece of drug data
nodrug<-read_excel("C:/model/nodrugexpt.xlsx")
## 1 dose, different ratios
Expt1 <-read_excel("C:/model/ratios.xlsx")
## 1 dose, different concs, fixed 50:50 ratio
Expt2 <-read_excel("C:/model/1dose.xlsx")
## 2 doses, different concs, fixed 50:50 ratio 
Expt3 <-read_excel("C:/model/2doses.xlsx")
## 3 doses, different concs, fixed 50:50 ratio
Expt4 <-read_excel("C:/model/3doses.xlsx")

# 2.Growth Model ####
growthmodel <- function(t,state,parms){
  with(as.list(c(state, parms)),{
    
    ##multiple doses
    pulse = c(rep(c(rep(1, mdur), rep(0, (24-mdur))),rounds),
              rep(0, (rounds*70*24)))
    
    delta_s = delta_r = 0;  
    #dose introduced within the first 3 days
    if(t<72){
      delta_s = approx(1:length(pulse),delta_s1*pulse,t, method="constant", rule=2)$y
      delta_r = approx(1:length(pulse),delta_r1*pulse,t, method="constant", rule=2)$y}
    
    ##convert PMR (parasite multiplication rates) to growth rates
    gs=(1/48)*log(Mr/(1-f))
    gr=(1/48)*log(Mr)
    
    ##model structure
    dS = gs*S - (1-x)*delta_s*S -x*delta_s*S+z*Sd   #growth of S
    dSd = x*delta_s*S-z*Sd                          #dormancy of S
    dR = gr*R - (1-y)*delta_r*R -y*delta_r*R+u*Rd   #growth of R
    dRd = y*delta_r*R-u*Rd                          #dormancy of R
    
    output <- c(dS,dSd,dR,dRd)
    list(output)
  })
}


# 3.Initial conditions and time ####
ratio<-0.58 #initial ratio in experiments

init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)

times <- seq(0, 70*24, by = 1) #in hours

# 4. Estimating fitness cost from no drug data####
nodrug_data<-nodrug$`% Wildype`
ratio<-0.5
init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)
#calculate least square of %C580Y
getssq<-function(parms){
  parms <- c(
    Mr=3.5,  #multiplication rate of R strain calculated from data
    f=parms[1],  #fitness cost of resistance
    delta_s1 = 0,#/h, drug killing effect on S strain
    delta_r1 = 0,#/h, drug killing effect on R strain
    x = 0 , #%dormancy in S
    y = 0, #%dormancy in R
    rounds = 1,     #rounds of doses
    mdur = 0.25*24,   #h,duration of drug treatment
    z=0.0000,
    u=0.0000
  )
  
  output_df <- data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
  total_parasite<-round(c(output_df$S+output_df$Sd+output_df$R+output_df$Rd ))
  perc_R <-(output_df$R+output_df$Rd)/total_parasite*100
  total_parasite<-total_parasite[seq(1,length(total_parasite),24)]
  perc_R <-perc_R[seq(24,length(perc_R),24)]
  t_nodrug<-c(1:11,14,17)
  perc_R<-perc_R[t_nodrug]
  
  ssq<-sum((perc_R-nodrug_data)^2)
  
  return(ssq)
  
}

fitssq<-optim(0.01,getssq,method="Brent",lower=0,upper=1)
f_fit<-fitssq$par

parms <- c(
  Mr=3.5,  #multiplication rate of R strain calculated from data
  f=f_fit,  #fitness cost of resistance
  delta_s1 = 0,#/h, drug killing effect on S strain
  delta_r1 = 0,#/h, drug killing effect on R strain
  x = 0 , #%dormancy in S
  y = 0, #%dormancy in R
  rounds = 1,     #rounds of doses
  mdur = 0.25*24,   #h,duration of drug treatment
  z=0.0000,
  u=0.0000
)
output_df <- data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
total_parasite<-round(c(output_df$S+output_df$Sd+output_df$R+output_df$Rd ))
perc_R <-(output_df$R+output_df$Rd)/total_parasite*100
total_parasite<-total_parasite[seq(1,length(total_parasite),24)]
perc_R <-perc_R[seq(24,length(perc_R),24)]
t_nodrug<-c(1:11,14,17)
perc_R<-perc_R[t_nodrug]
#plot(perc_R)
plot(t_nodrug,nodrug_data,type="p",pch=20, main="%C580Y in absence of drug using estimated f" ,xlab="Time(Days)",ylab="%C580Y")
lines(t_nodrug,perc_R,col="red")
print(f_fit)
# 5.Calculate NLL for growth model ####
r<-c(0.58,1/100,1/150,1/250,1/500,1/750,1/1000) #ratios
a<-0 #store sumNLL for each ratio 

getNLL<-function(parms){
  parms <- c(
    Mr=3.5,  #multiplication rate of R strain calculated from data
    f=0.2581026,  #fitness cost 
    delta_s1 = exp(parms[1]),#/h, drug clearance effect on S strain
    delta_r1 = exp(parms[2]),#/h, drug clearance effect on R strain
    x = exp(parms[3])/(1+exp(parms[3])), #% entering dormancy in S strain
    y = exp(parms[4])/(1+exp(parms[4])), #% entering dormancy in R strain
    rounds = 1,     #rounds of doses
    mdur =0.25*24,  #h,duration of drug treatment
    z=exp(parms[5]), #rate of waking up from dormancy in S strain
    u=exp(parms[6]) #rate of waking up from dormancy in R strain
  )
  
  for(i in 1:7){
    ratio<-r[i]
    if(i==1){
      expt_days<-Expt2$Days
    } else {
      expt_days<-Expt1$Days
    }
    init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)
    output_df <- data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
    total_parasite<-round(c(output_df$S+output_df$Sd+output_df$R+output_df$Rd ))
    perc_R <-(output_df$R+output_df$Rd)/total_parasite*100
    total_parasite<-total_parasite[seq(1,length(total_parasite),24)]
    perc_R <-perc_R[seq(1,length(perc_R),24)]
    
    total_parasite <-total_parasite[expt_days]
    perc_R <-perc_R[expt_days]

    if(i==1){
      expt_perc<-Expt2$`% C580Y...3`
      expt_tp<-Expt2$`Total Parasites...2`
      NLL<- -dpois(expt_perc,perc_R,log=TRUE)
      impossibles<- NLL== Inf
      NLL[impossibles]<- 999999
      NLL2<- -dpois(log10(expt_tp),log10(total_parasite),log=TRUE)
      impossibles<- NLL2== Inf
      NLL2[impossibles]<- 999999
      a[i]<-sum(NLL+4*NLL2)
    } else{
      expt_perc<-pull(Expt1[,i])
      NLL<- -dpois(expt_perc,perc_R,log=TRUE)
      impossibles<- NLL== Inf
      NLL[impossibles]<- 999999
      a[i]<-sum(NLL)
    }
    #print(expt_perc)
    
  }
  return(sumNLL<-sum(a))
}


# 6.Step1:Model fitting ####
fit1PoissonLL<-optim(c(log(3),log(2),log(0.15/(1-0.15)),log(0.2/(1-0.2)),log(0.0000002),log(0.002)),fn = getNLL,control=list(reltol=1e-12),hessian=TRUE )
sd<-sqrt(diag(solve(fit1PoissonLL$hessian)))
fitpar<-fit1PoissonLL$par
# sd<-sqrt(diag(solve(fit_model2$hessian)))
# fitpar<-fit_model2$par
fitparlq<-fitpar-1.96*sd #95% CI
fitparuq<-fitpar+1.96*sd #95% CI

fitpar<-cbind(fitpar,fitparlq,fitparuq)
fitpar[1,]<-exp(fitpar[1,])
fitpar[2,]<-exp(fitpar[2,])
fitpar[3,]<-exp(fitpar[3,])/(1+exp(fitpar[3,]))
fitpar[4,]<-exp(fitpar[4,])/(1+exp(fitpar[4,]))
fitpar[5,]<-exp(fitpar[5,])
fitpar[6,]<-exp(fitpar[6,])

print(fitpar)
# 7.Simulations with fitted parameters ####
fit<-NULL
for(i in 1:7){
    ratio<-r[i]
    if(i==1){
      expt_days<-Expt2$Days
      expt_perc<-Expt2$`% C580Y...3`
      expt_tp <-Expt2$`Total Parasites...2`
    } else {
      expt_days<-Expt1$Days
      expt_perc<-pull(Expt1[,i])
      expt_tp <-NULL
    }
    for (j in 1:3){
      parms <- c(
        Mr=3.5,  #multiplication rate of R strain
        f=0.2581026,  #fitness cost
        delta_s1 = fitpar[[1,j]],#/h, drug killing effect on S strain
        delta_r1 = fitpar[[2,j]],#/h, drug killing effect on R strain
        x = fitpar[[3,j]],
        y = fitpar[[4,j]],
        rounds = 1,     #rounds of doses
        mdur =0.25*24,  #h,duration of drug treatment
        z=fitpar[[5,j]],
        u=fitpar[[6,j]]
      )
    
    init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)
    
    fit.new<-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")
    fit<-cbind(fit,fit.new)
    fit.df<-data.frame(fit)
  #print(fit.df)
  #print(parms)
}}

# 8:Plot graphs ####
df1<-fit.df
dailydf1 <- df1[seq(1,nrow(df1),24),]
for(i in seq(1,105,by=5)){
  dailydf1[,i]<-(dailydf1[,i]-1)/24+1
}
t<- dailydf1$time

  ##Plot 1: 50:50
  i=1
  totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
  perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
  tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
  tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
  perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
  perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
  dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))

    Expt2<-as.data.frame(Expt2)
    ggplot(Expt2)+
      aes(x=Days,y=`Total Parasites...2`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=totalp1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=tp_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=tp_uq))+
      labs(title="Ratio 50:50 Total parasite vs Time",x="Days",y="Log Number of parasites")+
      scale_y_continuous(trans="log10")+
      theme_bw()
      
    ggplot(Expt2)+
      aes(x=Days ,y=`% C580Y...3`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 50:50: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))

    ##plot 2: 1:99
    i=16
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:99`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:99: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    ##plot 3: 1:149
    i=31
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:149`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:149: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    ##plot 4: 1:249
    i=46
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:249`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:249: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    ##plot 5: 1:499
    i=61
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:499`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:499: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    ##plot 6: 1:749
    i=76
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:749`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:749: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    ##plot 7: 1:999
    i=91
    totalp1 <- (dailydf1[,i+1] + dailydf1[,i+2] + dailydf1[,i+3]  + dailydf1[,i+4])
    perc_R1 <- (dailydf1[,i+3]  + dailydf1[,i+4])/totalp1*100
    tp_lq <-(dailydf1[,i+6] + dailydf1[,i+7] + dailydf1[,i+8]  + dailydf1[,i+9])
    tp_uq<-(dailydf1[,i+11] + dailydf1[,i+12] + dailydf1[,i+13]  + dailydf1[,i+14])
    perc_lq<-(dailydf1[,i+8]  + dailydf1[,i+9])/totalp1*100
    perc_uq<-(dailydf1[,i+13]  + dailydf1[,i+14])/totalp1*100
    dat<-data.frame(cbind(t,totalp1,perc_R1,tp_lq,perc_lq,tp_uq,perc_uq))
    
    Expt1<-as.data.frame(Expt1)
    ggplot(Expt1)+
      aes(x=Days,y=`1:999`)+
      geom_point()+
      geom_line(data=dat,colour="red",aes(x=t,y=perc_R1))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_lq))+
      geom_line(data=dat,colour="grey",aes(x=t,y=perc_uq))+
      labs(title="Ratio 1:999: %C580Y vs Time",x="Days",y="%C580Y")+
      theme_bw()+
      ylim(c(0,100))
    
    
print(fitpar)    
    
# 9:Step2:Second fitting####
fitpar2_store <-NULL
for (m in c(4,6,8,10)){
  ratio<-0.58
  init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)
  expt_days<-Expt2$Days
  expt_tp<-Expt2[[m]]
  expt_perc<-Expt2[[m+1]]
  na<-which(is.na(expt_perc))
  if(length(na)>0){
    expt_perc <-expt_perc[-na]
  }
  
  getNLL2<-function(parms){
    parms <- c(
      Mr=3.5,  #multiplication rate of R strain calculated from data
      f=0.2581026,  #fitness cost of resistance
      delta_s1 = exp(parms[1]),#/h, drug clearance effect on S strain
      delta_r1 = exp(parms[2]),#/h, drug clearance effect on R strain
      x = exp(parms[3])/(1+exp(parms[3])), #%dormancy in S
      y = exp(parms[4])/(1+exp(parms[4])), #%dormancy in R
      rounds = 1,     #rounds of doses
      mdur =0.25*24,  #h,duration of drug treatment
      z=exp(parms[5]),
      u=exp(parms[6])

    )
    output_df <- data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
    total_parasite<-round(c(output_df$S+output_df$Sd+output_df$R+output_df$Rd ))
    perc_R <-(output_df$R+output_df$Rd)/total_parasite*100
    total_parasite<-total_parasite[seq(1,length(total_parasite),24)]
    total_parasite <-total_parasite[expt_days]
    perc_R <-perc_R[seq(1,length(perc_R),24)]
    perc_R <-perc_R[expt_days]
    if(length(na)>0){
      perc_R <-perc_R[-na]
    }
    
    NLL<- -dpois(expt_perc,perc_R,log=TRUE)
    
    impossibles<- NLL== Inf
    NLL[impossibles]<- 999999
    sumNLL<-sum(NLL)
    
    NLL2<- -dpois(log10(expt_tp),log10(total_parasite), log=TRUE)
    
    impossibles2<- NLL2== Inf
    NLL2[impossibles2]<- 999999
    
    sumNLL<-sum(NLL+1000*NLL2)
    
    return(sumNLL)
  }
  
  # fitting outcome
  fit2PoissonLL<-optim(c(log(2.849260e+01/(2^((m-2)/2))),log(2.185359e+00/(2^((m-2)/2))),log(5.080380e-03/(1-5.080380e-03)),log(5.160803e-02/(1-5.160803e-02)),log(5.349312e-09),log(1.722116e-06)),fn = getNLL2,control=list(reltol=1e-12),hessian=TRUE )
  #fit2PoissonLL<-optim(c(log(2.74334940/(2^((m-2)/2))),log(1.53688635/(2^((m-2)/2))),log(0.04013134/(1-0.04013134)),log(0.42336408/(1-0.42336408 ))),fn = getNLL2,control=list(reltol=1e-12),hessian=TRUE )
  fitpar2<-fit2PoissonLL$par
  fitpar2[1]<-exp(fitpar2[1])
  fitpar2[2]<-exp(fitpar2[2])
  fitpar2[3]<-exp(fitpar2[3])/(1+exp(fitpar2[3]))
  fitpar2[4]<-exp(fitpar2[4])/(1+exp(fitpar2[4]))
  fitpar2[5]<-exp(fitpar2[5])
  fitpar2[6]<-exp(fitpar2[6])
  fitpar2_store<-cbind(fitpar2_store,fitpar2)
}
for(i in 1:4){
  parms <- c(
    Mr=3.5,  #multiplication rate of R strain
    f=0.2581026,  #fitness cost
    delta_s1 = fitpar2_store[[1,i]],#/h, drug killing effect on S strain
    delta_r1 = fitpar2_store[[2,i]],#/h, drug killing effect on R strain
    x = fitpar2_store[[3,i]],
    y = fitpar2_store[[4,i]],
    rounds = 1,     #rounds of doses
    mdur =0.25*24,   #h,duration of drug treatment
    z=fitpar2_store[[5,i]],
    u=fitpar2_store[[6,i]]
  )

  fit.df2<-data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))

  df2<-fit.df2
  dailydf2 <- df2[seq(1,nrow(df2),24),]
  dailydf2$time<-dailydf2$time/24+1
  t<- dailydf2$time
  totalp2 <- (dailydf2$S + dailydf2$Sd + dailydf2$R + dailydf2$Rd)
  perc_R2 <- (dailydf2$R + dailydf2$Rd)/totalp2*100
  dat2<-data.frame(cbind(t,totalp2,perc_R2)) 

  Expt2<-as.data.frame(Expt2)
  
  print(ggplot(Expt2)+
    aes(x=Days,y=Expt2[,i*2+2])+
    geom_point()+
    geom_line(data=dat2,colour="red",aes(x=t,y=totalp2))+
    labs(title="Ratio 50:50: Total parasite vs Time",x="Days",y="Log Number of parasites")+
    scale_y_continuous(trans="log10")+
    theme_bw())
  
  print(ggplot(Expt2)+
    aes(x=Days ,y=Expt2[,i*2+3])+
    geom_point()+
    geom_line(data=dat2,colour="red",aes(x=t,y=perc_R2))+
    labs(title="Ratio 50:50: %C580Y vs Time",x="Days",y="%C580Y")+
    theme_bw()+
    ylim(c(0,100)))
  print(m)
  print(fitpar2)
}


#### 10:Validation using Expt3,4####
parms <- c(
  Mr = 3.5,  #multiplication rate of R strain calculated from data
  f = 0.2581026,  #fitness cost
  delta_s1 = 2.849260e+01,#/h, drug clearance effect on S strain
  delta_r1 = 2.185359e+00,#/h, drug clearance effect on R strain
  x = 5.080380e-03, #% entering dormancy in S strain
  y = 5.160803e-02, #% entering dormancy in R strain
  rounds = 1,     #rounds of doses
  mdur = 0.25*24,  #h,duration of drug treatment
  z = 5.349312e-09, #rate of waking up from dormancy in S strain
  u = 1.722116e-06 #rate of waking up from dormancy in R strain
)

init_state <- c(S = 1.32*10^7, Sd=0, R = 1.68*10^7, Rd=0)
# 100%DHA ####
# Model1:1 dose, 100% DHA, 50:50 ratio
model1 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode") #vode uses a variable step-size for calculation, so fast change can be captured
model<-model1
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=2
Expt<-as.data.frame(Expt2)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="1 dose, 100%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="1 dose, 100%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  

# Model12: 2 doses, 100%DHA, 50:50

parms[7]<-2

model12 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model12

df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=2
Expt<-as.data.frame(Expt3)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="2 doses, 100%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="2 doses, 100%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  

# Model13: 3 doses, 100%DHA, 50:50

parms[7]<-3

model13 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model13
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=2
Expt<-as.data.frame(Expt4)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="3 doses, 100%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="3 doses, 100%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  



# 1 dose, different ratios ####
# Model3: 1 dose, 100% DHA, 1:99 ratio 
parms[7]<-1
init_state <- c(S = 3*10^7*0.99, Sd=0, R = 3*10^7*0.01, Rd=0)

model3 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")
model<-model3
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=2
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  

# Model7: 1 dose, 100% DHA, 1:149 ratio

init_state <- c(S = 3*10^7*(149/150), Sd=0, R = 3*10^7*(1/150), Rd=0)

model7 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model7
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=3
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  

# Model8: 1 dose, 100% DHA, 1:249 ratio 

init_state <- c(S = 3*10^7*(249/250), Sd=0, R = 3*10^7*(1/250), Rd=0)

model8 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")
model<-model8
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=4
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  

# Model9: 1 dose, 100% DHA, 1:499 ratio 

init_state <- c(S = 3*10^7*(499/500), Sd=0, R = 3*10^7*(1/500), Rd=0)

model9 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")
model<-model9
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=5
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  


# Model10: 1 dose, 100% DHA, 1:749 ratio 

init_state <- c(S = 3*10^7*(749/750), Sd=0, R = 3*10^7*(1/750), Rd=0)

model10 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model10
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=6
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))  


# Model11: 1 dose, 100% DHA, 1:999 ratio 

init_state <- c(S = 3*10^7*(999/1000), Sd=0, R = 3*10^7*(1/1000), Rd=0)

model11 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model11
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=7
Expt<-as.data.frame(Expt1)
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="%C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100)) 


# 50%DHA ####
# Model2:1 dose,50% DHA,50:50 ratios

parms[3]<-3.36E+01
parms[4]<-2.38E-01
parms[5]<-1.66E-01
parms[6]<-1.07E-02
parms[9]<-1.30E-09
parms[10]<-5.48E-06
init_state <- c(S = 1.32*10^7, Sd=0, R = 1.68*10^7, Rd=0)
model2 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model2
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=4
Expt<-as.data.frame(Expt2)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="1 dose, 50%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="1 dose, 50%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# Model14: 2 doses, 50%DHA, 50:50

parms[7]<-2

model14 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model14
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=4
Expt<-as.data.frame(Expt3)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="2 doses, 50%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="2 doses, 50%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# Model15: 3 doses, 50%DHA, 50:50

parms[7]<-3

model15 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model15
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=4
Expt<-as.data.frame(Expt4)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="3 doses, 50%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="3 doses, 50%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))



# 25%DHA ####
# Model4:1 dose,25% DHA,50:50 ratios
parms[7]<-1
parms[3]<-3.56E+00
parms[4]<-2.73E-01
parms[5]<-5.08E-03
parms[6]<-5.16E-02
parms[9]<-5.35E-09
parms[10]<-1.72E-06
init_state <- c(S = 1.32*10^7, Sd=0, R = 1.68*10^7, Rd=0)
model4 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model4
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=6
Expt<-as.data.frame(Expt2)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="1 dose, 25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="1 dose, 25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# Model16: 2 doses, 25%DHA, 50:50

parms[7]<-2

model16 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model16
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=6
Expt<-as.data.frame(Expt3)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="2 doses, 25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="2 doses, 25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))
# Model17: 3 doses, 25%DHA, 50:50

parms[7]<-3

model17 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model17
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=6
Expt<-as.data.frame(Expt4)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="3 doses, 25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="3 doses, 25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# 12.5%DHA ####
# Model5:1 dose,12.5% DHA,50:50 ratios
parms[7]<-1
parms[3]<-1.23E+00
parms[4]<-8.71E-02
parms[5]<-1.81E-02
parms[6]<-5.81E-02
parms[9]<-1.40E-08
parms[10]<-1.15E-06
init_state <- c(S = 1.32*10^7, Sd=0, R = 1.68*10^7, Rd=0)
model5 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model5
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=8
Expt<-as.data.frame(Expt2)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="1 dose, 12.5%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="1 dose, 12.5%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# Model18: 2 doses, 12.5%DHA, 50:50

parms[7]<-2

model18 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model18
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=8
Expt<-as.data.frame(Expt3)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="2 doses, 12.5%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="2 doses, 12.5%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))


# Model19: 3 doses, 12.5%DHA, 50:50

parms[7]<-3

model19 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model19
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=8
Expt<-as.data.frame(Expt4)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="3 doses, 12.5%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="3 doses, 12.5%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# 6.25%DHA ####
# Model6:1 dose,6.25% DHA,50:50 ratios
parms[7]<-1
parms[3]<-1.06E+00
parms[4]<-4.51E-01
parms[5]<-2.46E-03
parms[6]<-4.47E-02
parms[9]<-9.90E-09
parms[10]<-1.96E-06
init_state <- c(S = 1.32*10^7, Sd=0, R = 1.68*10^7, Rd=0)
model6 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model6
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=10
Expt<-as.data.frame(Expt2)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="1 dose,6.25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="1 dose,6.25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))

# Model20: 2 doses, 6.25%DHA, 50:50

parms[7]<-2

model20 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model20
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=10
Expt<-as.data.frame(Expt3)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="2 doses,6.25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="2 doses,6.25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))


# Model21: 3 doses, 6.25%DHA, 50:50

parms[7]<-3

model21 <-ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode")

model<-model21
df3 <- as.data.frame(model)
dailydf3 <- df3[seq(1,nrow(df3),24),]
dailydf3$time<-dailydf3$time/24
t3<- dailydf3$time+1
totalp3 <- (dailydf3$S + dailydf3$Sd + dailydf3$R + dailydf3$Rd)
perc_R3 <- (dailydf3$R + dailydf3$Rd)/totalp3*100
dat3 <-as.data.frame(cbind(t3,totalp3,perc_R3))
m=10
Expt<-as.data.frame(Expt4)
#Total parasite
ggplot(Expt)+
  aes(x=Days,y=Expt[,m])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=totalp3))+
  labs(title="3 doses,6.25%DHA : Total parasite vs Time",x="Days",y="Log Number of parasites")+
  scale_y_continuous(trans="log10")+
  theme_bw()
#%C580Y
ggplot(Expt)+
  aes(x=Days ,y=Expt[,m+1])+
  geom_point()+
  geom_line(data=dat3,colour="red",aes(x=t,y=perc_R3))+
  labs(title="3 doses,6.25%DHA : %C580Y vs Time",x="Days",y="%C580Y")+
  theme_bw()+
  ylim(c(0,100))



#11:Sensitivity Analysis####
#parameter values from first fitting
parms <- c(
  Mr=3.5,  #multiplication rate of R strain calculated from data
  f=0.2581026,  #fitness cost of resistance
  delta_s1 = 2.849260e+01, #/h, drug killing effect on S strain
  delta_r1 = 2.185359e+00,#/h, drug killing effect on R strain
  x = 5.080380e-03, #%dormancy in S
  y = 5.160803e-02, #%dormancy in R
  rounds = 1 ,  #rounds of doses
  mdur = 0.25*24 ,
  z= 5.349312e-09,
  u= 1.722116e-06
)

times <- seq(1, 70*24, by = 1)
ratio<-0.58
init_state <- c(S = 3*10^7*(1-ratio), Sd=0, R = 3*10^7*ratio, Rd=0)
output_df <- data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
total_parasite<-round(c(output_df$S+output_df$Sd+output_df$R+output_df$Rd ))
perc_R <-(output_df$R+output_df$Rd)/total_parasite*100
total_parasite<-total_parasite[seq(1,length(total_parasite),24)]
perc_R <-perc_R[seq(1,length(perc_R),24)]
plot(Expt2$Days,Expt2$`% C580Y...3`,type="p",pch=20)
lines(perc_R,type="l",col="red")
output_df

x<-c(0.005*2,0.005*0.5,0.005*10,0.005/10)
y<-c(0.05*2,0.05*0.5,0.05*10,0.05/10)
z<-c((5*10^-9)*2,(5*10^-9)*0.5,(5*10^-9)*10,(5*10^-9)/10)
u<-c((2*10^-6)*2,(2*10^-6)*0.5,(2*10^-6)*10,(2*10^-6)/10)
# x<-c(0.001,0.005,0.01,0.05,0.1)
# y<-c(0.001,0.01,0.05,0.1,0.5)
# z<-c(0.000000001,0.000000005,0.00000001,0.0000001,0.00001)
# u<-c(0.0000001,0.000001,0.00001,0.0001,0.001)
par(mar=c(5, 4, 4, 8))
plot(Expt2$Days,Expt2$`% C580Y...3`,type="p",pch=20,main="Sensitivity analysis for x",xlab="Days",ylab="%C580Y")
lines(perc_R,type="l",col="red")
for(i in 1:4){
  colour<-c("green","blue","pink","orange")
  parms[5]<-x[i]
  sa_xout<-data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
  sa_xtp<-round(c(sa_xout$S+sa_xout$Sd+sa_xout$R+sa_xout$Rd ))
  sa_xperc<-(sa_xout$R+sa_xout$Rd)/sa_xtp*100
  sa_xtp<-sa_xtp[seq(1,length(sa_xtp),24)]
  sa_xperc <-sa_xperc[seq(1,length(sa_xperc),24)]
  # dat4<-as.data.frame(cbind(t,sa_xtp,sa_xperc))
  # sa_plot_new<-sa_plot_new+geom_line(data=dat4,aes(x=t,y=sa_xperc))
  # sa_plot_new
  lines(sa_xperc,type="l",col=colour[i])
  legend("topright", inset=c(-0.3, 0), legend=c("10","2","0.5","0.1","baseline"),
         col=c("pink","green","blue","orange","red"), lty=1,cex=0.8,xpd=TRUE)
}

plot(Expt2$Days,Expt2$`% C580Y...3`,type="p",pch=20,main="Sensitivity analysis for y",xlab="Days",ylab="%C580Y")
lines(perc_R,type="l",col="red")
for(i in 1:4){
  colour<-c("green","blue","pink","orange")
  parms[5]<-5.080380e-03
  parms[6]<-y[i]
  sa_xout<-data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
  sa_xtp<-round(c(sa_xout$S+sa_xout$Sd+sa_xout$R+sa_xout$Rd ))
  sa_xperc<-(sa_xout$R+sa_xout$Rd)/sa_xtp*100
  sa_xtp<-sa_xtp[seq(1,length(sa_xtp),24)]
  sa_xperc <-sa_xperc[seq(1,length(sa_xperc),24)]
  lines(sa_xperc,type="l",col=colour[i])
  legend("topright", inset=c(-0.3, 0), legend=c("10","2","0.5","0.1","baseline"),
         col=c("pink","green","blue","orange","red"), lty=1,cex=0.8,xpd=TRUE)
}

plot(Expt2$Days,Expt2$`% C580Y...3`,type="p",pch=20, main="Sensitivity analysis for z",xlab="Days",ylab="%C580Y")
lines(perc_R,type="l",col="red")
for(i in 1:4){
  colour<-c("green","blue","pink","orange")
  parms[5]<-5.080380e-03
  parms[6]<-5.160803e-02
  parms[9]<-z[i]
  sa_xout<-data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
  sa_xtp<-round(c(sa_xout$S+sa_xout$Sd+sa_xout$R+sa_xout$Rd ))
  sa_xperc<-(sa_xout$R+sa_xout$Rd)/sa_xtp*100
  sa_xtp<-sa_xtp[seq(1,length(sa_xtp),24)]
  sa_xperc <-sa_xperc[seq(1,length(sa_xperc),24)]
  lines(sa_xperc,type="l",col=colour[i])
  legend("topright", inset=c(-0.3, 0), legend=c("10","2","0.5","0.1","baseline"),
         col=c("pink","green","blue","orange","red"), lty=1,cex=0.8,xpd=TRUE)
}
plot(Expt2$Days,Expt2$`% C580Y...3`,type="p",pch=20,main="Sensitivity analysis for u",xlab="Days",ylab="%C580Y")
lines(perc_R,type="l",col="red")
for(i in 1:4){
  colour<-c("green","blue","pink","orange")
  parms[5]<-5.080380e-03
  parms[6]<-5.160803e-02
  parms[9]<-5.349312e-09
  parms[10]<-u[i]
  sa_xout<-data.frame(ode(times=times, y=init_state, func=growthmodel,parms=parms,method="vode"))
  sa_xtp<-round(c(sa_xout$S+sa_xout$Sd+sa_xout$R+sa_xout$Rd ))
  sa_xperc<-(sa_xout$R+sa_xout$Rd)/sa_xtp*100
  sa_xtp<-sa_xtp[seq(1,length(sa_xtp),24)]
  sa_xperc <-sa_xperc[seq(1,length(sa_xperc),24)]
  lines(sa_xperc,type="l",col=colour[i])
  legend("topright", inset=c(-0.3, 0), legend=c("10","2","0.5","0.1","baseline"),
         col=c("pink","green","blue","orange","red"), lty=1,cex=0.8,xpd=TRUE)
}
