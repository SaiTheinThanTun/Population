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

#reading in files for census data
setwd("C:/wd")
source('avg_stk_tbl.R')
age_prob_0to97 <- read.csv("C:/wd/0to97_age_prob.csv", header=FALSE)
male_prob_0to97 <- read.csv("C:/wd/0to97_male_prob.csv", header=FALSE)
#age, 98+ were unaccounted for
age <- 0:97

####parameters####
durinf <- 7 #duration of infection ###may need to readjust when transforming into shiny
a <- .1 #human blood feeding rate
b <- .3 #probability of disease transmission per bite for human
c <- .7 #probability a mosquito becomes infected after biting an infected human
muo <- .05 ##10 days survival= 20 half-days survival, therefore 1/20=.05
mui <- .05

H <- 80 #human population
X <- 10 #infected humans
M <- 800 #initial mosquito population
Z <- 200 #initial infected mosquitos
timesteps_days <- 1095 #28
timeres <- 1 #time resolution of .5 days
timesteps <- timesteps_days*(1/timeres) #365*2 
no_sims <- 10 #50 #10 #no. of simulations

lci <- .05 #lowest confidence interval
hci <- .95  #highest confidence interval


recover <- timeres/durinf #1/(2*durinf) #probability of getting recovered


####synthesizing age and gender####
sim_age <- sample(age, H, replace=TRUE,prob=age_prob_0to97[,1])
#hist(sim_age) 
gender <- rep(NA,length(sim_age))

for(i in 1:length(sim_age)){
  p <- male_prob_0to97[sim_age[i]+1,1] #male_prob is already arranged in ascending age
  gender[i] <- sample(2,1,prob=c(p,1-p))
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

no.patch.x <- 4 #no. of patches across x
no.patch.y <- 4 #no. of patches across y
total.patch <- no.patch.x*no.patch.y
random_no <- random_no2 <- patch <- rep(NA, H)


df <- as.data.frame(cbind(sim_age,gender,infected_h, random_no, random_no2, patch)) #variable addition for populated dataframe

####codebook for df, modify this after finalizing####
#1. sim_age
#2. gender
#3. infected_h
#4. tts #time to become susceptible again
#5. random_n #random no. drawn from uniform distribution
#6. random_n2 # ...
#7. current #a switch to detect if an individual is infected in current timestep or not


#######outputting csv: write an initialized file#####
write.csv(df, file='0.csv')


#######Simulate Summary table function#####
simulate_summ <- function(){#function for subsequent timesteps
  
  summ_tab <- matrix(NA, nrow=timesteps+1, ncol=8) # summary table for plotting, +1 because it starts from 0 #variable addition for simulation table
  colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam', 'lam_h2') #column names for the summary table
  #variable addition for simulation table

  summ_tab[,1] <- seq(0,timesteps_days,by=timeres)
  
  for(j in 0:timesteps){ #this means 2:(timesteps+1)
    seas <- (sin(.01722*timeres*j)*.02)+.1
    
      
    #this also needs to be changed during the subsequent timesteps
    m <- M/H
    z <- Z/M
    X <- sum(df$infected_h) #no. of infected humans
    S <- (M-Z) #susceptible mosquitos
    x <- X/H #ratio of infectious humans
    #rate of change of Z from ODE
    lam <- a*c*x #1-(1-(a*c))^x #a*c*x ###Reed-Frost
    lam_h2 <- m*a*b*z
    
    
    
    
    
    for(i in 1:nrow(df)){
      
      
      df$random_no[i] <- runif(1) #drawing random no. for each individual
      df$random_no2[i] <- runif(1)
      df$patch[i] <- sample(total.patch,1)
      lam_h <- seas*(sum(df[which(df$patch==df$patch[i]),]$infected_h)/length(which(df$patch==df$patch[i]))) #infectedpersonsSamePatch

      if(df$infected_h[i]==0){
        if(df$random_no[i]<=lam_h){
          df$infected_h[i] <- 1
        }
      } else{
        if(df$random_no2[i]<=recover){
          df$infected_h[i] <- 0
        }
      }
    }
    #writing a summary table
    #summ_tab[j,1] <- j
    k <- j+1 #because the loop starts from 0
    summ_tab[k,2] <- H-X
    summ_tab[k,3] <- X
    summ_tab[k,4] <-lam_h
    summ_tab[k,5] <- S #############################
    summ_tab[k,6] <- Z #need to have some limitation on Z, infected mosquitos
    summ_tab[k,7] <- lam
    summ_tab[k,8] <- lam_h2
    
    S_prev <- S
    
    S <- S_prev+M*mui-muo*S_prev-lam*S_prev
    Z <- Z+lam*S_prev-muo*Z
    
    M <- S+Z #recalculating mosquito population
    
    ######outputting csv of the simulation on each timestep#######
    if(j<20 | j>(max(timesteps)-10)){
      write.csv(df, file=paste(j,".csv",sep=""))
    }
  }
  summ_tab
}

summ_tab <- simulate_summ() #this is to be used for plotting a single simulation


####plotting 1 simulation####
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

####plotting multiple simulation####


#creating a list of simulations
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
