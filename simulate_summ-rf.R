#######Simulate Summary table function#####
#lambda_h adapted with Reed-Frost
simulate_summ <- function(){#function for subsequent timesteps
  
  summ_tab <- matrix(NA, nrow=timesteps+1, ncol=7) # summary table for plotting, +1 because it starts from 0 #variable addition for simulation table
  colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h','S','Z','lam') #column names for the summary table
  #variable addition for simulation table
  summ_tab[1,] <- time0 #the first line of the table. the states at time0
  
  #there's an error which one to take as time 0 (or 0.5)
  summ_tab[,1] <- seq(0,timesteps_days,by=(1/2))
  
  for(j in 1:timesteps+1){ #this means 2:(timesteps+1)
    
    for(i in 1:nrow(df)){
      if(df[i,5]<=lam_h){ #if uniform random no. drawn for individual is <= prob of infected
        df[i,3] <- df[i,6] <- 1 #denoting this person is infected on this timestep
      }
      
      if(df[i,3]==1 && df[i,6]==1){ #if infected #at current timestep 
        
        df[i,4] <- rnorm(1,mean=1,sd=.2) * durinf #input into tts, time to become susceptable again
        
      }
      
      df[i,4] <- df[i,4]-.5 #tts-.5 per timestep
      
      if(df[i,4]<=0 && df[i,3]==1){ #currently infected, but durinf is over
        df[i,3] <- 0 #then he becomes suscepitable again
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
    lam <- 1-(1-(a*c))^x #a*c*x
    Z <- Z+lam*(M-Z)
    #m <- M/H ###no. of mosquitos doesn't change FOR NOW
    z <- Z/M
    lam_h <- 1-(1-(a*b))^(m*z) #m*a*b*z
    
    #writing a summary table
    #summ_tab[j,1] <- j
    summ_tab[j,2] <- H-X
    summ_tab[j,3] <- X
    summ_tab[j,4] <-lam_h
    summ_tab[j,5] <- M-Z 
    summ_tab[j,6] <- Z #need to have some limitation on Z, infected mosquitos
    summ_tab[j,7] <- lam
    
    ######outputing csv of the simulation on each timestep#######
    #if(j<10 | j>(max(timesteps)-10)){
    #  write.csv(df, file=paste(j,".csv",sep=""))
    #}
  }
  summ_tab
}