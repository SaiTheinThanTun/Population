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
timesteps_days <- 28
timesteps <- timesteps_days*2 #365*2 

#this also needs to be changed during the subsequent timesteps
m <- M/H
z <- Z/M
x <- X/H #ratio of infectious humans

lam_h <- m*a*b*z
#lam_h <- 0.075

lam <- a*c*x #lambda for mosquitos


#synthesizing a population
#pop_size <- 5000

#age, 98+ were unaccounted for
age <- 0:97



sim_age <- sample(age, H, replace=TRUE,prob=age_prob_0to97[,1])
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
infected_h <- rep(NA,H)
for(i in 1:H){
  infected_h[i] <- sample(c(0,1),1, prob=c(.8,.2))
}

tts <- rep(0, H) #time to become susceptable, 1/dur_inf in normal distribution

random_no <- runif(H) # random uniform no. to decide the prob. of being infected if susceptable
current <- rep(1, H) #infected in current timestep

df <- cbind(sim_age,gender,infected_h, tts, random_no, current)

#codebook for df
#1. sim_age
#2. gender
#3. infected_h
#4. tts #time to become susceptible again
#5. random_n #random no. drawn from uniform distribution
#6. current #a switch to detect if an individual is infected in current timestep or not


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
summ_tab <- matrix(NA, nrow=timesteps, ncol=7) # summary table for plotting
colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h','M','Z','lam')

#there's an error which one to take as time 0 (or 0.5)
summ_tab[,1] <- seq(0.5,timesteps_days,by=(1/2))
  
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
  X <- sum(df[,3]) #no. of infected humans
  x <- X/H #ratio of infectious humans
  #rate of change of Z from ODE
  lam <- a*c*x
  Z <- Z+lam*(M-Z)
  
  #m <- M/H ###no. of mosquitos doesn't change FOR NOW
  
  z <- Z/M
  lam_h <- m*a*b*z
  
  #writing a summary table
  #summ_tab[j,1] <- j
  summ_tab[j,2] <- H-X
  summ_tab[j,3] <- X
  summ_tab[j,4] <-lam_h
  summ_tab[j,5] <- M 
  summ_tab[j,6] <- Z #need to have some limitation on Z, infected mosquitos
  summ_tab[j,7] <- lam
  
  if(j<10 | j>(max(timesteps)-10)){
    write.csv(df, file=paste(j,".csv",sep=""))
  }
}


#plotting
par(mar=c(5,4,4,4))
plot(summ_tab[,1],summ_tab[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop with lambda ",lam_h()))
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible humans",side=2,line=2.5) 

box()
par(new=TRUE)
plot(summ_tab[,1],summ_tab[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
axis(4, ylim=c(0,17),col="red") 
mtext("Infected humans",side=4, line=2.5)

axis(1,pretty(range(summ_tab[,1]),10))
mtext("Time (0.5 days)",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))


#infected var
#1 <= infected
#0 <- susceptible

#synthesizing patches of mosquitos

write.csv(summ_tab,file=paste('summary_ibm_',Sys.Date(),'.csv',sep=''))
          