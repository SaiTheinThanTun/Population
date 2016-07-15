#just a snippet of code to compare speed

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
a <- .15 #human blood feeding rate
b <- .3 #probability of disease transmission per bite for human
c_ <- .7 #probability a mosquito becomes infected after biting an infected human
muo <- .05 ##10 days survival= 20 half-days survival, therefore 1/20=.05
mui <- .05

seas_switch <- 1 #logical switch for seasonality
amp <- .2
phi <- 0
magnitude <- 1


H <- 300 #human population
X <- 20 #infected humans
M <- 800 #800 #initial mosquito population
Z <- 200 #200 #initial infected mosquitos
timesteps_days <- 365 #730 # 1095 #28
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

no.patch.x <- 1 #no. of patches across x
no.patch.y <- 1 #no. of patches across y
total.patch <- no.patch.x*no.patch.y

prob_infected <- prob_recovery <- patch.lam_m <- patch.lam_h <- patch <- random_no <- random_no2 <- rep(NA, H)

#patch.mosq <- matrix(NA,total.patch,3) #init mosq data on each patch
patch.mosq <- matrix(rep(c(Z,M, M-Z),total.patch),total.patch,3, byrow = TRUE)

df <- as.data.frame(cbind(sim_age,gender,infected_h, random_no, random_no2, patch, patch.lam_h, patch.lam_m,prob_infected, prob_recovery)) #variable addition for populated dataframe

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
  colnames(summ_tab) <- c('timesteps','susceptables_avg','infected_avg', 'seas','S_avg','Z_avg','lam_m_vector_avg', 'lam_h_vector_avg') #column names for the summary table
  #variable addition for simulation table
  
  summ_tab[,1] <- seq(0,timesteps_days,by=timeres)
  
  for(j in 0:timesteps){ #this means 2:(timesteps+1)
    seas <- amp*cos(2*pi*((j*timeres)-phi)/365)+magnitude #(sin(.01722*timeres*j)*.02)+.2
    seas <- (seas*seas_switch)+(1-seas_switch)
    #(value*switch) + (1-value)
    
    df$patch <- sample(total.patch,H, replace=TRUE)
    df$patch <- as.factor(df$patch) #,levels=as.character(1:total.patch))
    if(length(levels(df$patch))!=total.patch) print("some patch(es) have no human")
    df$infected_h <- as.factor(df$infected_h)
    df$random_no <- runif(H)
    df$random_no2 <- runif(H)
    
    
    patch.h <- as.matrix(table(df$patch,df$infected_h))
    if(length(levels(df$infected_h))==1){
      if(levels(df$infected_h)==0){ #all uninfected
        X <- 0
        H_patch <- patch.h[,1]
      }
      else { #all infected
        X <- H_patch <- patch.h[,1]
      }
    }
    else { #normal procedure
      X <- patch.h[,2] 
      H_patch <- (patch.h[,1]+X)
    }
    #     X <- patch.h[,2] 
    #     H_patch <- (patch.h[,1]+X)
    x <- X/H_patch
    
    #this also needs to be changed during the subsequent timesteps
    Z <- patch.mosq[as.numeric(levels(df$patch)),1]
    M <- patch.mosq[as.numeric(levels(df$patch)),2]
    m <- M/H_patch # M/H, vector for m
    z <- Z/M #Z/M, vector for z
    S <- patch.mosq[as.numeric(levels(df$patch)),3] #M-Z #susceptible mosquitos
    
    #X <- sum(df$infected_h) #no. of infected humans
    #S <- (M-Z) #susceptible mosquitos
    #x <- X/H #ratio of infectious humans
    #rate of change of Z from ODE
    
    
    #as df$patch is a vector, there;s no need to check of the consistency with the total.patch
    ####resolving NA and NaN values
    x[is.na(x)] <- 0
    x[is.nan(x)] <- 0
    m[is.na(m)] <- 0
    m[is.nan(m)] <- 0
    z[is.na(z)] <- 0
    z[is.nan(z)] <- 0
    #either above or check the validity of lam_h and lam_m vectors themselves
    
    #     x <- as.vector(by(df$infected_h,df$patch,sum) / by(df$infected_h,df$patch,length)) #this is x by patch
    #     if(!all(1:total.patch %in% unique(df$patch))){ #this solves the situation where no individual is on a particular patch
    #       putback0 <- which(!(1:total.patch %in% df$patch))
    #       
    #       for(k in 0:(length(putback0)-1)){
    #         x <- append(x,0,after=putback0[k+1]-1)
    #       }
    #     } 
    
    lam_m_vector <- a*c_*x*seas  #1-(1-(a*c))^x #a*c*x ###Reed-Frost seasonality into lam_m
    lam_h_vector <- m*a*b*z
    #seas*x #seas*(sum(df[which(df$patch==df$patch[i]),]$infected_h)/length(which(df$patch==df$patch[i])))
    H_summ <- length(df$infected_h)
    X_summ <- sum(df$infected_h==1)
    #writing a summary table ###maybe move this up
    #summ_tab[j,1] <- j
    k <- j+1 #because the loop starts from 0
    summ_tab[k,2] <- H_summ-X_summ #median(H_patch-X)
    summ_tab[k,3] <- X_summ #median(X)
    summ_tab[k,4] <- seas
    summ_tab[k,5] <- mean(S) #############################
    summ_tab[k,6] <- mean(Z) #need to have some limitation on Z, infected mosquitos
    summ_tab[k,7] <- mean(lam_m_vector)
    summ_tab[k,8] <- mean(df$prob_infected) #, na.rm=TRUE) #mean(prob_recovery) #mean(lam_h_vector)
    
    
    
    for(i in 1:nrow(df)){
      
      #lam_h <- seas*(sum(df[which(df$patch==df$patch[i]),]$infected_h)/length(which(df$patch==df$patch[i]))) #infectedpersonsSamePatch
      df$patch.lam_h[i] <- lam_h_vector[df$patch[i]]
      df$patch.lam_m[i] <- lam_m_vector[df$patch[i]]
      
      df$prob_infected[i] <- 1-exp(-df$patch.lam_h[i]*timeres)
      df$prob_recovery[i] <- 1-exp(-recover*timeres)
      
      if(df$infected_h[i]==0){ #if not infected
        if(df$random_no[i]<= df$prob_infected[i]){ #if getting infected #1-exp(-k*t)
          #this calculated for each patch and use here efficiently
          df$infected_h[i] <- 1
        }
      } else{
        if(df$random_no2[i]<= df$prob_recovery[i]){ #if infected, but getting recovered
          df$infected_h[i] <- 0
        }
      }
    }
    
    
    S_prev <- NA
    #for loop for patches
    for(i2 in 1:total.patch){
      S_prev[i2] <- patch.mosq[i2,3]
      
      
      patch.mosq[i2,3] <- S_prev[i2]+M[i2]*mui-muo*S_prev[i2]-lam_m_vector[i2]*S_prev[i2]
      patch.mosq[i2,1] <- Z[i2]+lam_m_vector[i2]*S_prev[i2]-muo*Z[i2]
      patch.mosq[i2,2] <- patch.mosq[i2,1] + patch.mosq[i2,3]
      #       S_prev[i2] <- S[i2]
      #       
      #       S[i2] <- S_prev[i2]+M[i2]*mui-muo*S_prev[i2]-lam_m_vector[i2]*S_prev[i2]
      #       Z[i2] <- Z[i2]+lam_m_vector[i2]*S_prev[i2]-muo*Z[i2]
      #       
      #       M[i2] <- S[i2]+Z[i2]
    }
    #recalculating mosquito population
    
    ######outputting csv of the simulation on each timestep#######
    if(k<20 | k>(max(timesteps)-10)){
      write.csv(df, file=paste(k,".csv",sep=""))
    }
  }
  summ_tab
}

#summ_tab <- simulate_summ() #this is to be used for plotting a single simulation
system.time( simulate_summ() )
