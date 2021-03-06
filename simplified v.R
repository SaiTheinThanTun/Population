#before you run, create a directory named 'wd' under 'C:/'
#copy a csv file: china_age_prob.csv into the 'C:/wd'


#reading in files for census data
setwd("C:/wd")
age_prob <- read.csv("C:/wd/china_age_prob.csv", header=FALSE)[,1]
age <- 0:(length(age_prob)-1)

####parameters####
durinf <- 7 #duration of infection ###may need to readjust when transforming into shiny

amp <- .2
phi <- 210
magnitude <- .8


H <- 80 #human population
X <- 40 #infected humans


timesteps <- 1095 #365*2 



recover <- 1/durinf #1/(2*durinf) #probability of getting recovered


####synthesizing age and gender####
sim_age <- sample(age, H, replace=TRUE,prob=age_prob)
#hist(sim_age) 
#gender <- rep(NA,length(sim_age))
prob_male <- .55

gender <- sample(c(0,1),H, replace=TRUE,prob=c(1-prob_male, prob_male))


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
patch <- sample(total.patch,H, replace=TRUE)
random_no <- random_no2 <- rep(NA, H)


df <- as.data.frame(cbind(sim_age,gender,infected_h, random_no, random_no2, patch)) #variable addition for populated dataframe




#######Simulate Summary table function#####
simulate_summ <- function(){#function for subsequent timesteps
  
  summ_tab <- matrix(NA, nrow=timesteps+1, ncol=4) # summary table for plotting, +1 because it starts from 0 #variable addition for simulation table
  colnames(summ_tab) <- c('timesteps','susceptables','infected', 'lam_h') #column names for the summary table
  
  #summ_tab[,1] <- seq(0,timesteps,by=1)
  
  for(j in 0:timesteps){ #this means 2:(timesteps+1)
    seas <- amp*cos(2*pi*(j-phi)/365)+magnitude #(sin(.01722*j)*.02)+.2
    
    
    #this also needs to be changed during the subsequent timesteps
#     m <- M/H
#     z <- Z/M
     X <- sum(df$infected_h) #no. of infected humans
#     S <- (M-Z) #susceptible mosquitos
#     x <- X/H #ratio of infectious humans
#     #rate of change of Z from ODE
#     lam <- a*c*x #1-(1-(a*c))^x #a*c*x ###Reed-Frost
#     lam_h2 <- m*a*b*z
    
    preval.patch <- as.vector(by(df$infected_h,df$patch,sum) / by(df$infected_h,df$patch,length))
    if(!all(1:total.patch %in% unique(df$patch))){ #this solves the situation where no individual is on a particular patch
      putback0 <- which(!(1:total.patch %in% df$patch))
      
      for(k in 0:(length(putback0)-1)){
        preval.patch <- append(preval.patch,0,after=putback0[k+1]-1)
      }
    } 
    
    lam_h_vector <- seas*preval.patch #seas*(sum(df[which(df$patch==df$patch[i]),]$infected_h)/length(which(df$patch==df$patch[i])))
    
    
    
    
    for(i in 1:nrow(df)){
      
      
      df$random_no[i] <- runif(1) #drawing random no. for each individual
      df$random_no2[i] <- runif(1)
      df$patch[i] <- sample(total.patch,1)
      #lam_h <- seas*(sum(df[which(df$patch==df$patch[i]),]$infected_h)/length(which(df$patch==df$patch[i]))) #infectedpersonsSamePatch
      lam_h <- lam_h_vector[df$patch[i]]
      
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
    k <- j+1 #because the loop starts from 0
    summ_tab[k,1] <- j
    summ_tab[k,2] <- H-X
    summ_tab[k,3] <- X
    summ_tab[k,4] <-lam_h
    
    
    
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
mtext(paste("Timesteps (Time resolution: ",1," day)",sep=""),side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

###writing csv 1 simulation####
write.csv(summ_tab,file=paste('summary_ibm_',Sys.Date(),'.csv',sep=''))


