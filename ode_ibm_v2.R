#before you run, create a directory named 'wd' under 'C:/'
#copy 2 csv files: 0to97_age_prob.csv and 0to97_male_prob.csv into the 'C:/wd'
#copy a generic function file called 'avg_stk_tbl.R' into 'C:/wd/'

#this R script file does the following:
#1. synthesize a population based upon Myanmar census data 2014, see ####codebook for df#### section for variables
#2. input parameters are required to synthesize and simulate the population, see ####parameters#### section
#3. create a function that runs a single simulation and summarize the result on each timestep
#4. run multiple simulations. you can change #no. of simulations in ####parameters####
####and create their averages with confidence interval
#5. write csv files and plot the results

#if more variables are to be added, search for this hashtag #variable addition
#if output is needed for each timestep, search for and comment out this hashtag #outputting csv

####libraries####
#library(ggplot2)
library(deSolve)

#reading in files for census data
setwd("C:/wd")
source('avg_stk_tbl.R')
age_prob <- read.csv("C:/wd/0to97_age_prob.csv", header=FALSE)[,1]
male_prob_0to97 <- read.csv("C:/wd/0to97_male_prob.csv", header=FALSE) #consider this later
#age, 98+ were unaccounted for
age <- 0:(length(age_prob)-1)

####parameters####
durinf <- 7 #duration of infection ###may need to readjust when transforming into shiny
a <- .5 #human blood feeding rate
b <- .3 #probability of disease transmission per bite for human
c_ <- .7 #probability a mosquito becomes infected after biting an infected human
muo <- 0 #.05 ##10 days survival= 20 half-days survival, therefore 1/20=.05
mui <- 0 #.05

seas_switch <- 1 #logical switch for seasonality
amp <- .2
phi <- 0
magnitude <- 1 #.8


H <- 300 #human population
X <- 1 #infected humans
M <- 80 #800 #initial mosquito population
Z <- 20 #200 #initial infected mosquitos
timesteps_days <- 365 #1095 #28
timeres <- 1 #time resolution of .5 days
timesteps <- timesteps_days*(1/timeres) #365*2 
no_sims <- 10 #50 #10 #no. of simulations

lci <- .05 #lowest confidence interval
hci <- .95  #highest confidence interval


recover <- timeres/durinf #1/(2*durinf) #probability of getting recovered


####synthesizing age and gender####
sim_age <- sample(age, H, replace=TRUE,prob=age_prob)
#hist(sim_age) 
gender <- rep(NA,length(sim_age))

for(i in 1:length(sim_age)){
  p <- male_prob_0to97[sim_age[i]+1,1] #male_prob is already arranged in ascending age
  gender[i] <- sample(c(0,1),1,prob=c(1-p,p))
}

####testing the proportions####
#tmp <- cbind(sim_age,gender)
#tmp <- as.data.frame(tmp)
#colnames(tmp) <- c('s_age','s_gender')
#tmp$s_gender <- as.factor(tmp$s_gender)
#qplot(s_age, data=tmp,fill=s_gender)

####Initializing / synthesizing infected population####
infected_h <- rep(NA,H)
S_prob <- (H-X)/H
I_prob <- X/H
for(i in 1:H){
  infected_h[i] <- sample(c(0,1),1, prob=c(S_prob,I_prob))
}

no.patch.x <- 1 #4 #no. of patches across x
no.patch.y <- 1 #4 #no. of patches across y
total.patch <- no.patch.x*no.patch.y

prob_infected <- prob_recovery <- patch.lam_m <- patch.lam_h <- patch <- random_no <- random_no2 <- rep(NA, H)

#patch.mosq <- matrix(NA,total.patch,3) #init mosq data on each patch
patch.mosq <- matrix(rep(c(Z,M, M-Z),total.patch),total.patch,3, byrow = TRUE)

df <- as.data.frame(cbind(sim_age,gender,infected_h, random_no, random_no2, patch, patch.lam_h, patch.lam_m,prob_infected, prob_recovery)) #variable addition for populated dataframe


summ_tab <- simulate_summ()

par(mar=c(5,4,4,4))
plot(summ_tab[,1],summ_tab[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop")) # with lambda",lam_h))
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible humans",side=2,line=2.5) 

box()
par(new=TRUE)
plot(summ_tab[,1],summ_tab[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
axis(4, ylim=c(0,17),col="red") 
mtext("Infected humans",side=4, line=2.5)

axis(1,pretty(range(summ_tab[,1]),10))
mtext(paste("Timesteps (Time resolution: ",timeres," day)",sep=""),side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

###writing csv 1 simulation####
write.csv(summ_tab,file=paste('summary_ibm_',Sys.Date(),'.csv',sep=''))

####multiple simulation####
sims <- list()
for(i in 1:no_sims){
  sims[[i]] <- simulate_summ()
}

#averaging across the list
avg_sims <- avg_stk_tbl(sims)

# #lower CI (LCI)
lci_sims <- avg_stk_tbl(sims,'ci',ci=lci)

# #high CI (HCI)
hci_sims <- avg_stk_tbl(sims,'ci',ci=hci)
colnames(avg_sims) <- colnames(hci_sims) <- colnames(lci_sims) <- colnames(summ_tab) #c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table

par(mar=c(5,4,4,4))
plot(avg_sims[,1],avg_sims[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop")) # with lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,2], rev(lci_sims[,2])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible humans",side=2,line=2.5) 

box()
par(new=TRUE)
plot(avg_sims[,1],avg_sims[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,3], rev(lci_sims[,3])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
axis(4, ylim=c(0,17),col="red") 
mtext("Infected humans",side=4, line=2.5)

axis(1,pretty(range(avg_sims[,1]),10))
mtext("Time",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

#mosquitos
par(mar=c(5,4,4,4))
plot(avg_sims[,1],avg_sims[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("mosquito_pop")) # with lambda",lam_h,"and CI",lci,'-',hci))
polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,5], rev(lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible mosquitos",side=2,line=2.5) 

box()
par(new=TRUE)
plot(avg_sims[,1],avg_sims[,6], type="l", col="red", axes=FALSE, xlab="", ylab="")
polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,6], rev(lci_sims[,6])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
axis(4, ylim=c(0,17),col="red") 
mtext("Infected mosquitos",side=4, line=2.5)

axis(1,pretty(range(avg_sims[,1]),10))
mtext("Time",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

###writing csv_ average of multiple simulations####
write.csv(avg_sims,file=paste('avg_summary_ibm_',Sys.Date(),'.csv',sep=''))


###end of IBM####

parameters <- c(
  c(mui=mui,    # birth #lifespan of mosquito 10 days
    muo=muo,    # death
    #beta= #per capita effective contact with infected human per unit time
    #ce = (.3*.01), #probability of disease transmission per bite * biting rate
    a = a, #human blood feeding rate
    b = b, #probability of disease transmission per bite for human
    c = c_, ##probability of disease transmission per bite for mosquitos
    #beta = input$beta, #probability of disease transmission per bite for mosquitos
    recover = recover, #duration of infection in days * 2, because of 1/2 day time step
    #immunity = input$immunity * 2 #duration of immunity * 2, because of 1/2 day time step
    seas_switch = seas_switch, #logical switch for seasonality
    amp = amp,
    phi = phi,
    magnitude = magnitude
  ))


  lam_h <- (M/H)*a*b*(Z/M)



  initSh <- H-X  #+initRh) #no of suscepitables
  
  initS<-M-Z
  initD <- 0
  
  ##lam_h <- 0
  lam <- 0
  #seas_0 <- ((amp*cos(2*pi*(0-phi)/365)+magnitude)*seas_switch)+(1-seas_switch)
  
  state <- c(Sh= initSh, X = X, lam_h = lam_h, S = initS, Z = Z, lam=lam,Y=0)
  times <- seq(0,timesteps, by=timeres)
  
  # set up a function to solve the model
  mosQ<-function(t, state, parameters) 
  {
    with(as.list(c(state, parameters)),
         {
           
           # define variables
           M <- (S+Z)
           H <- (Sh+X)
           m <- M/H #ratio of mosquitos to humans
           z <- Z/M #ratio of infectious mosquitos
           x <- X/H #ratio of infectious humans
           
           seas <- ((amp*cos(2*pi*(Y-phi)/365)+magnitude)*seas_switch)+(1-seas_switch) #(sin(.01722*timeres*j)*.02)+.2
           #(value*switch) + (1-value)
           #seas<-1+amp*cos(2*pi*(Y-phi)/52)
           #beta<-R0*(muo+nui)*gamma/(muo+gamma)
           lam <- a*c*x*seas
           lam_h <- m*a*b*z
           
           # rate of change for mosquitos
           dS <- mui*M-muo*S-lam*S
           dZ <- -muo*Z+lam*S
           #dD <- muo*S+muo*Z #remove this!!!!
           #dM <- 0
           
           # rate of change for humans
           dSh <- -lam_h*Sh+recover*X #durinf*2 because input was in days and ODE was in .5 days
           dX <- lam_h*Sh-recover*X  #durinf*2 because input was in days and ODE was in .5 days
           #dRh <- (1/durinf)*X-(1/immunity)*Rh
           dY <- 1
           
           # return the rate of change
           list(c(dSh, dX,lam_h, dS, dZ, lam, dY))
         }
    ) 
    
  }
  out <- ode(y = state, times = times, func = mosQ, parms = parameters)
  

    write.csv(out,paste('summary_ode_',Sys.Date(),'.csv',sep=''))


# 
#   par(mar=c(5,4,4,4)) #default is par(mar=c(5,4,4,2))
#   
  plot(out[,1],out[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="mosq_pop")
  axis(2, ylim=c(0,17),col="blue") 
  mtext("Susceptible mosquitoes",side=2,line=2.5) 
  box()
  par(new=TRUE)
  plot(out[,1],out[,6], type="l", col="red", axes=FALSE, xlab="", ylab="")
  axis(4, ylim=c(0,17),col="red") 
  mtext("Infected mosquitoes",side=4,line=2.5)
  
  axis(1,pretty(range(out[,1]),10))
  mtext("Time (days)",side=1,col="black",line=2.5)
  
  legend("top",legend=c("Susceptibles","Infected"),
         text.col=c("blue","red"),pch= "__", col=c("blue","red"))

    #humans ode
    par(mar=c(5,4,4,4))
    plot(out[,1],out[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="human_pop with lambda ")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible humans",side=2,line=2.5) 
    
    box()
    par(new=TRUE)
    plot(out[,1],out[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected humans",side=4, line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (days)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
    
    
    ###comparison graph####
    #plot from avg_sims and out
    
    par(mar=c(5,4,4,4))
    plot(avg_sims[,1],avg_sims[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop")) # with lambda",lam_h,"and CI",lci,'-',hci))
    lines(out[,1],out[,2], lty = 2, col="darkblue", lwd=3)
    polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,2], rev(lci_sims[,2])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible humans",side=2,line=2.5) 
    
    box()
    par(new=TRUE)
    plot(avg_sims[,1],avg_sims[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
    lines(out[,1],out[,3], lty = 2, col="darkred", lwd=3)
    polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,3], rev(lci_sims[,3])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected humans",side=4, line=2.5)
    
    axis(1,pretty(range(avg_sims[,1]),10))
    mtext("Time",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
    
    #mosquitos
    par(mar=c(5,4,4,4))
    plot(avg_sims[,1],avg_sims[,5], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("mosquito_pop")) # with lambda",lam_h,"and CI",lci,'-',hci))
    lines(out[,1],out[,5], lty = 2, col="darkblue", lwd=3)
    polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,5], rev(lci_sims[,5])),col=rgb(0,0,100,50,maxColorValue=255), border=NA)
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible mosquitos",side=2,line=2.5) 
    
    box()
    par(new=TRUE)
    plot(avg_sims[,1],avg_sims[,6], type="l", col="red", axes=FALSE, xlab="", ylab="")
    lines(out[,1],out[,6], lty = 2, col="darkred", lwd=3)
    polygon(c(avg_sims[,1], rev(avg_sims[,1])), c(hci_sims[,6], rev(lci_sims[,6])),col=rgb(100,0,0,50,maxColorValue = 255), border=NA)
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected mosquitos",side=4, line=2.5)
    
    axis(1,pretty(range(avg_sims[,1]),10))
    mtext("Time",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))