#before you run, create a directory named 'wd' under 'C:/'
#copy 2 csv files: 0to97_age_prob.csv and 0to97_male_prob.csv into the 'C:/wd'

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
library(ggplot2)

#reading in files for census data
setwd("C:/wd")
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
X <- 30 #infected humans
M <- 800 #initial mosquito population
Z <- 200 #initial infected mosquitos
timesteps_days <- 28
timesteps <- timesteps_days*2 #365*2 
no_sims <- 10 #no. of simulations

lci <- .05 #lowest confidence interval
hci <- .95  #highest confidence interval

#this also needs to be changed during the subsequent timesteps
m <- M/H
z <- Z/M
x <- X/H #ratio of infectious humans

lam_h <- m*a*b*z #lam_h <- 0.075 #lambda for humans
lam <- a*c*x #lambda for mosquitos


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

####synthesizing infected population####
infected_h <- rep(NA,H)
S_prob <- (H-X)/H
I_prob <- X/H
for(i in 1:H){
  infected_h[i] <- sample(c(0,1),1, prob=c(S_prob,I_prob))
}

tts <- rep(0, H) #time to become susceptable, 1/dur_inf in normal distribution

random_no <- runif(H) # random uniform no. to decide the prob. of being infected if susceptable
current <- rep(1, H) #infected in current timestep

df <- cbind(sim_age,gender,infected_h, tts, random_no, current) #variable addition for populated dataframe

####codebook for df####
#1. sim_age
#2. gender
#3. infected_h
#4. tts #time to become susceptible again
#5. random_n #random no. drawn from uniform distribution
#6. current #a switch to detect if an individual is infected in current timestep or not


###initializing####
for(i in 1:nrow(df)){
  if(df[i,3] && df[i,6]){ #if infected #at current timestep 
    
    
    df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to susceptable
    
    
  }
  
  df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
  df[i,6] <- 0 # resetting 'infected at current timestep'
}

#first row of the summary table
#time 0
X <- sum(df[,3]) #no. of infected humans
x <- X/H #ratio of infectious humans
#rate of change of Z from ODE
lam <- a*c*x #1-(1-(a*c))^x #a*c*x ###Reed-Frost
S_prev <- (M-Z)
Z_prev <- Z #only in this first step, on later steps in the function, it was not recorded
M_prev <- M
S <- S_prev+M_prev*mui-muo*S_prev-lam*S_prev
Z <- Z_prev+lam*S_prev-muo*Z_prev

M <- S+Z #recalculating mosquito population
###need to check this####

#m <- M/H ###no. of mosquitos doesn't change FOR NOW
z <- Z/M
lam_h <- m*a*b*z

time0 <- c(0, H-X, X, lam_h, S_prev, Z_prev, lam) #variable addition for simulation table
##above, lam_h and lam values are for the next time step

#######outputting csv: write an initialized file#####
write.csv(df, file='0.csv')


#######Simulate Summary table function#####
simulate_summ <- function(){#function for subsequent timesteps
  
  summ_tab <- matrix(NA, nrow=timesteps+1, ncol=7) # summary table for plotting, +1 because it starts from 0 #variable addition for simulation table
  colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table
  #variable addition for simulation table
  summ_tab[1,] <- time0 #the first line of the table. the states at time0
  
  #there's an error which one to take as time 0 (or 0.5)
  summ_tab[,1] <- seq(0,timesteps_days,by=(1/2))
  
  for(j in 1:timesteps+1){ #this means 2:(timesteps+1)
    
    for(i in 1:nrow(df)){
      if(df[i,3]==0 & df[i,5]<=lam_h){ #if uniform random no. drawn for 'uninfected' individual is <= prob of getting infected
        df[i,3] <- df[i,6] <- 1 #denoting this person is infected on this timestep
      }
      
      if(df[i,3]==1 && df[i,6]==1){ #if infected #at current timestep 
        
        df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to become susceptable again
        
      }
      
      df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
      
      if(df[i,4]<=0 && df[i,3]==1){ #currently infected, but durinf is over
        df[i,3] <- 0 #then he becomes suscepitable again on the next timestep
      }
      
      #resetting for the next round
      df[i,5] <- runif(1) #drawing random no. for each individual
      df[i,6] <- 0 # resetting 'infected at current timestep'
    }
    #at the end of big for loop
    #calculate summary variables and lam_h for the next timestep
    X <- sum(df[,3]) #no. of infected humans
    x <- X/H #ratio of infectious humans
    #rate of change of Z from ODE
    lam <- a*c*x #1-(1-(a*c))^x #a*c*x ###Reed-Frost
    S_prev <- (M-Z)
    S <- S_prev+M*mui-muo*S_prev-lam*S_prev
    Z <- Z+lam*S_prev-muo*Z
    
    M <- S+Z #recalculating mosquito population
    #m <- M/H ###no. of mosquitos doesn't change FOR NOW
    z <- Z/M
    lam_h <- m*a*b*z #1-(1-(a*b*m))^z #m*a*b*z  ###Reed-Frost
    
    #writing a summary table
    #summ_tab[j,1] <- j
    summ_tab[j,2] <- H-X
    summ_tab[j,3] <- X
    summ_tab[j,4] <-lam_h
    summ_tab[j,5] <- S #############################
    summ_tab[j,6] <- Z #need to have some limitation on Z, infected mosquitos
    summ_tab[j,7] <- lam
    
    ######outputting csv of the simulation on each timestep#######
    if(j<20 | j>(max(timesteps)-10)){
      write.csv(df, file=paste(j-1,".csv",sep=""))
    }
  }
  summ_tab
}

summ_tab <- simulate_summ() #this is to be used for plotting a single simulation


####plotting 1 simulation####
par(mar=c(5,4,4,4))
plot(summ_tab[,1],summ_tab[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop with lambda",lam_h))
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

###writing csv 1 simulation####
write.csv(summ_tab,file=paste('summary_ibm_',Sys.Date(),'.csv',sep=''))

####plotting multiple simulation####


#creating a list of simulations
sims <- list()
for(i in 1:no_sims){
  sims[[i]] <- simulate_summ()
}

#averaging across the list
tmp_avg <- rep(NA,no_sims)
avg_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_avg[k] <- sims[[k]][j,i]
    }
    avg_sims[j,i] <- mean(tmp_avg)
  }
}

#lower CI (LCI)
tmp_lci <- rep(NA,no_sims)
lci_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_lci[k] <- sims[[k]][j,i]
    }
    lci_sims[j,i] <- quantile(tmp_lci, probs=lci, na.rm=TRUE)
  }
}

#high CI (HCI)
tmp_hci <- rep(NA,no_sims)
hci_sims <- matrix(NA,nrow(sims[[1]]),ncol(sims[[1]])) #initializing a blank dataset of summary table
for(i in 1:ncol(sims[[1]])){ #outer loop for the columns
  for(j in 1:nrow(sims[[1]])){ #inner loop for the rows
    for(k in 1:no_sims){#innermost loop for no. of simulations(3rd dimension)
      tmp_hci[k] <- sims[[k]][j,i]
    }
    hci_sims[j,i] <- quantile(tmp_hci, probs = hci, na.rm=TRUE)
  }
}

colnames(avg_sims) <- colnames(hci_sims) <- colnames(lci_sims) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table

par(mar=c(5,4,4,4))
plot(avg_sims[,1],avg_sims[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("human_pop with lambda",lam_h,"and CI",lci,'-',hci))
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
mtext("Time (0.5 days)",side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

###writing csv_ average of multiple simulations####
write.csv(avg_sims,file=paste('avg_summary_ibm_',Sys.Date(),'.csv',sep=''))
