#synthesizing a population
pop_size <- 5000

#age, 98+ were unaccounted for
age <- 0:97

age_prob_0to97 <- read.csv("D:/Dropbox/IBM project_Sai/Population/0to97_age_prob.csv", header=FALSE)
#age_prob_0to97 <- read.csv("C:/Users/lisa/Dropbox/IBM project_Sai/Population/0to97_age_prob.csv", header=FALSE)


sim_age <- sample(age, pop_size, replace=TRUE,prob=age_prob_0to97[,1])
hist(sim_age) 

#gender
gender <- rep(NA,length(sim_age))

male_prob_0to97 <- read.csv("D:/Dropbox/IBM project_Sai/Population/0to97_male_prob.csv", header=FALSE)
#male_prob_0to97 <- read.csv("C:/Users/lisa/Dropbox/IBM project_Sai/Population/0to97_male_prob.csv", header=FALSE)
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


#subsequent timesteps
for(i in 1:nrow(df)){
  if(df[i,5]>lam_h){ #if uniform random no. drawn for individual is > prob of infected
    df[i,3] <- df[i,6] <- 1 #denoting this person is infected on this timestep
  }
  
  if(df[i,3] && df[i,6]){ #if infected #at current timestep 
    
    
    df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to susceptable
    
    
  }
  
  df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
  df[i,5] <- runif(1) #drawing random no. for each individual
  df[i,6] <- 0 # resetting 'infected at current timestep'
}
#at the end of big for loop
#calculate lam_h

if('S'){
  if(runif() > prob.ds)
  {
    thatman <- 'X'
  }
} else
{
  
}


#infected var
#1 <= infected
#0 <- susceptible

#synthesizing patches of mosquitos

