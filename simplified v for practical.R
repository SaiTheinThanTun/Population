#before you run, create a directory named 'wd' under 'C:/'
#copy a csv file: age_prob.csv into the 'C:/wd'

library(ggplot2)
library(plyr)

setwd("C:/wd")
age_prob <- read.csv("C:/wd/age_prob.csv", header=FALSE)[,1] #reading in file for age distribution

####parameters####
P <- 200 #total population
I <- 190 #initial infected individuals

durinf <- 7 #duration of infection ###may need to readjust when transforming into shiny
recover <- 1/durinf #1/(2*durinf) #probability of getting recovered
amp <- .2
phi <- 210 #peak of the seasonal wave (eg at 210th day of 365)
magnitude <- .4 #
timesteps <- 365 #1095 #365*2 

prob_male <- .55 #probability of being male
no.patch.x <- 4 #no. of patches across x
no.patch.y <- 4 #no. of patches across y


####synthesizing a population####
age <- sample(0:(length(age_prob)-1), P, replace=TRUE,prob=age_prob) #create a vector of age 0 to the last age in the age_prob file, and sample from it with age_prob values
male <- sample(c(0,1),P, replace=TRUE,prob=c(1-prob_male, prob_male))
infected <- c(rep(1,I),rep(0,P-I))
total.patch <- no.patch.x*no.patch.y
patch.lam <- random_no <- random_no2 <- patch <- rep(NA, P)

df <- as.data.frame(cbind(age,male,infected, random_no, random_no2, patch, patch.lam))


#####plotting age pyramid####
df$male <- as.factor(df$male)
ggplot(data=df, aes(x=age,fill=male, colour=male))+
  geom_histogram(subset=.(male==1), binwidth=5, alpha=.8)+
  geom_histogram(subset=.(male==0),aes(y=..count..*(-1)),binwidth=5, alpha=.8)+
  scale_y_continuous(breaks=seq(-(dim(df)[1]*.2),dim(df)[1]*.2,5),labels=abs(seq(-(dim(df)[1]*.2),dim(df)[1]*.2,5)))+
  coord_flip()

##creating a table to store summary information for each timestep####
summ_tab <- matrix(NA, nrow=timesteps+1, ncol=3) # summary table for plotting, +1 because timestep will start from 0
colnames(summ_tab) <- c('timesteps','susceptables','infected') #column names for the summary table

for(j in 0:timesteps){ #this means 2:(timesteps+1)
  seas <- amp*cos(2*pi*(j-phi)/365)+magnitude 
  df$patch <- sample(total.patch,P, replace=TRUE)
  df$random_no <- runif(P)
  df$random_no2 <- runif(P)
  
  I <- sum(df$infected) #no. of infected
  
  preval.patch <- as.vector(by(df$infected,df$patch,sum) / by(df$infected,df$patch,length)) #calculating prevalance on each patch
  if(!all(1:total.patch %in% unique(df$patch))){ #this solves the situation where no individual is on a particular patch
    putback0 <- which(!(1:total.patch %in% df$patch))
    for(k in 0:(length(putback0)-1)){
      preval.patch <- append(preval.patch,0,after=putback0[k+1]-1)
    }
  } 
  
  lam_vector <- seas*preval.patch 
  
  for(i in 1:nrow(df)){
    df$patch.lam[i] <- lam_vector[df$patch[i]]
    
    if(df$infected[i]==0){
      if(df$random_no[i]<=df$patch.lam[i]){
        df$infected[i] <- 1
      }
    } else{
      if(df$random_no2[i]<=recover){
        df$infected[i] <- 0
      }
    }
  }
  
  #writing into the summary table
  k <- j+1 #because the loop starts from 0
  summ_tab[k,1] <- j
  summ_tab[k,2] <- P-I
  summ_tab[k,3] <- I
  
  ######output csv of the simulation on each timestep#######
  #This is gonna get your working directory a bit messy. But if you're curious about what's happening 
  #during the first and the last 10 timesteps, remove the # from the next line and re-run the whole scripts.
  #if(j<10 | j>(max(timesteps)-10)){    write.csv(df, file=paste(j,"timestep.csv",sep=""))    }
}

summ_tab

####plotting the simulation####
par(mar=c(5,4,4,4))
plot(summ_tab[,1],summ_tab[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main=paste("Plot of summary table"))
axis(2, ylim=c(0,17),col="blue") 
mtext("Susceptible",side=2,line=2.5) 

box()
par(new=TRUE)
plot(summ_tab[,1],summ_tab[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
axis(4, ylim=c(0,17),col="red") 
mtext("Infected",side=4, line=2.5)

axis(1,pretty(range(summ_tab[,1]),10))
mtext(paste("Timesteps (Time resolution: ",1," day)",sep=""),side=1,col="black",line=2.5)

legend("top",legend=c("Susceptibles","Infected"),
       text.col=c("blue","red"),pch= "__", col=c("blue","red"))

###writing csv file for simulation####
write.csv(summ_tab,file=paste('summary_ibm_',Sys.Date(),'.csv',sep=''))
