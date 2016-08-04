#ball game in IBM
initS <- 19 #initial susceptible population
initI1 <- 1 
beta <- 2 #.1 #effective contact
no.of.timesteps <- 10

#S -> I1 -> I2 -> R
#0 -> 1 -> 2 -> 3
pop <- append(rep(0,initS),rep(1,initI1))

#create table to store data
sim.table <- as.data.frame(matrix(NA, no.of.timesteps, 6))
colnames(sim.table) <- c('WeekNo','S','I1','I2','R','Total')
tmp <- list()
lambda <- NA
real <- NA
  
for(j in 1:no.of.timesteps){
#   lambda[j] <- 1-(1-beta)^sum(pop==1 | pop==2) #/length(pop)) #calculating probability of getting infected ##Reed-Frost
#   lambda[j] <- beta*sum(pop==1| pop==2) #/length(pop)
  lambda[j] <- 1-exp(-beta*sum(pop==1| pop==2)/length(pop))
  real[j] <- beta*sum(pop==1 | pop ==2)/length(pop)
  random.no <- runif(length(pop)) #creating random no.
  
  sim.table$WeekNo[j] <- j
  sim.table$S[j] <- sum(pop==0)
  sim.table$I1[j] <- sum(pop==1)
  sim.table$I2[j] <- sum(pop==2)
  sim.table$R[j] <- sum(pop==3)
  sim.table$Total[j] <- sim.table$S[j]+sim.table$I1[j]+sim.table$I2[j]+sim.table$R[j]
  tmp[[j]] <- pop #remove
  
  for(i in 1:length(pop)){
    if(pop[i]==0){
      if(random.no[i]<lambda[j]){
        pop[i] <- pop[i] + 1
      }
    }
    
    else if(pop[i]==1){
      pop[i] <- pop[i] + 1
    }
    
    else if(pop[i]==2){
      pop[i] <- pop[i] + 1
    }
  }
}

plot(sim.table$WeekNo, sim.table$S, type='l', col='blue', main='Ballgame in IBM', xlab='Weeks', ylab='No. of individuals', lwd=2)
lines(sim.table$WeekNo, sim.table$I1+sim.table$I2, type='l', col='red', lwd=2)
lines(sim.table$WeekNo, sim.table$R, type='l', lwd=2)
legend(sim.table$WeekNo[1],length(pop)/2,  c("S","I","R"), lty=c(1,1,1), lwd=c(2,2,2),col=c("blue","red","black"))

sim.table

# tmp <- do.call(rbind, tmp) #remove
# 
# write.csv(tmp, 'tmp.csv') #remove
#in terms of generic functions

# library(RCurl)
# 
# script <- getURL("https://github.com/SaiTheinThanTun/generic-functions/blob/9acfdd500019835c5e08f4680330589e8ea90385/syn_pop.R", ssl.verifypeer = FALSE)
# source(script)
# eval(parse(text = script))
source("D:\\Dropbox\\IBM project_Sai\\generic functions\\syn_pop.R")
source("D:\\Dropbox\\IBM project_Sai\\generic functions\\state_trans.R")
pop2 <- syn_pop(c(19,1,0,0))

beta <- 2 #.1 #effective contact
no.of.timesteps <- 10

#define lambda outside for loop for the first value
lambda = beta*sum(pop2[,2],pop2[,3])/sum(pop2)
sim.table <- matrix(NA, no.of.timesteps, ncol(pop2))
for(i in 1:10){ #for.ibm(eval)
   
  pop2 <- state.trans(3,4,100,pop2) 
  pop2 <- state.trans(2,3,100,pop2)
  pop2 <- state.trans(1,2,lambda,pop2) 
  
  lambda = beta*sum(pop2[,2],pop2[,3])/sum(pop2)
  sim.table[i,] <- colSums(pop2)
}
pop2
sim.table <- as.data.frame(cbind(1:no.of.timesteps,sim.table))
colnames(sim.table) <- c('WeekNo','S','I1','I2','R')
sim.table

plot(sim.table$WeekNo, sim.table$S, type='l', col='blue', main='Ballgame in IBM', xlab='Weeks', ylab='No. of individuals', lwd=2)
lines(sim.table$WeekNo, sim.table$I1+sim.table$I2, type='l', col='red', lwd=2)
lines(sim.table$WeekNo, sim.table$R, type='l', lwd=2)
legend(sim.table$WeekNo[1],length(pop)/2,  c("S","I","R"), lty=c(1,1,1), lwd=c(2,2,2),col=c("blue","red","black"))
