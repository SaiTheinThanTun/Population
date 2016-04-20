#to rearrange the codes into a structure
setwd("C:/wd")

#reading in files
age_prob_0to97 <- read.csv("D:/Dropbox/IBM project_Sai/Population/0to97_age_prob.csv", header=FALSE)
#age_prob_0to97 <- read.csv("C:/Users/lisa/Dropbox/IBM project_Sai/Population/0to97_age_prob.csv", header=FALSE)

male_prob_0to97 <- read.csv("D:/Dropbox/IBM project_Sai/Population/0to97_male_prob.csv", header=FALSE)
#male_prob_0to97 <- read.csv("C:/Users/lisa/Dropbox/IBM project_Sai/Population/0to97_male_prob.csv", header=FALSE)

#inputs
durinf <- 7 #duration of infection
a <- .1 #human blood feeding rate
b <- .3 #probability of disease transmission per bite for human
c <- .7 #probability a mosquito becomes infected after biting an infected human

H <- 80 #human population
X <- 30 #infected humans
M <- 800 #initial mosquito population
Z <- 200 #initial infected mosquitos
timesteps <- 28*2 #365*2 

#this also needs to be changed during the subsequent timesteps
m <- M/H
z <- Z/M
x <- X/H #ratio of infectious humans

lam_h <- m*a*b*z
#lam_h <- 0.075

lam <- a*c*x #lambda for mosquitos


#synthesizing a population
pop_size <- 5000

#age, 98+ were unaccounted for
age <- 0:97



sim_age <- sample(age, pop_size, replace=TRUE,prob=age_prob_0to97[,1])
hist(sim_age) 

#gender
gender <- rep(NA,length(sim_age))

for(i in 1:length(sim_age)){
  p <- male_prob_0to97[sim_age[i]+1,1] #male_prob is already arranged in ascending age
  gender[i] <- sample(2,1,prob=c(p,1-p))
}

#testing the proportions
library(ggplot2)
tmp <- cbind(sim_age,gender)
tmp <- as.data.frame(tmp)
colnames(tmp) <- c('s_age','s_gender')
tmp$s_gender <- as.factor(tmp$s_gender)
qplot(s_age, data=tmp,fill=s_gender)

#disease state
#malaria <- rep(NA,pop_size)
#for(i in 1:pop_size){
#  malaria[i] <- sample(c("S","I","R"),1, prob=c(.7,.1,.2))
#}
infected_h <- rep(NA,pop_size)
for(i in 1:pop_size){
  infected_h[i] <- sample(c(0,1),1, prob=c(.8,.2))
}

tts <- rep(0, pop_size) #time to become susceptable, 1/dur_inf in normal distribution

random_no <- runif(5000) # random uniform no. to decide the prob. of being infected if susceptable
current <- rep(1, pop_size) #infected in current timestep

df <- cbind(sim_age,gender,infected_h, tts, random_no, current)

#initializing
for(i in 1:nrow(df)){
  if(df[i,3] && df[i,6]){ #if infected #at current timestep 
    
    
    df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to susceptable
    
    
  }
  
  df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
  df[i,6] <- 0 # resetting 'infected at current timestep'
}
write.csv(df, file='0.csv')



#subsequent timesteps
summ_tab <-  # summary table for plotting
  
for(j in 1:timesteps){
  
  for(i in 1:nrow(df)){
    if(df[i,5]<=lam_h){ #if uniform random no. drawn for individual is <= prob of infected
      df[i,3] <- df[i,6] <- 1 #denoting this person is infected on this timestep
    }
    
    if(df[i,3]==1 && df[i,6]==1){ #if infected #at current timestep 
      
      df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to susceptable
      
    }
    
    df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
    
    if(df[i,4]<=0 && df[i,3]==1){ #infected at current step, but durinf is over
      df[i,3] <- 0
    }
    
    #resetting for the next round
    df[i,5] <- runif(1) #drawing random no. for each individual
    df[i,6] <- 0 # resetting 'infected at current timestep'
  }
  #at the end of big for loop
  #calculate lam_h
  X <- sum(df[,3])
  x <- X/H #ratio of infectious humans
  #rate of change of Z from ODE
  lam <- a*c*x
  Z <- Z+lam*(M-Z)
  
  #m <- M/H ###no. of mosquitos doesn't change FOR NOW
  
  z <- Z/M
  lam_h <- m*a*b*z
  
  if(j<10 | j>(max(timesteps)-10)){
    write.csv(df, file=paste(j,".csv",sep=""))
  }
}



#infected var
#1 <= infected
#0 <- susceptible

#synthesizing patches of mosquitos

