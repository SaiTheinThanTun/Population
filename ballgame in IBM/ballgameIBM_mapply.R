#Apply in the IBM####
initS <- 19
initI1 <- 1
beta <- 1.5 #effective contact
no.of.timesteps <- 10
pop <- append(rep(0,initS),rep(1,initI1))
rno <- runif(length(pop))

sim.table <- as.data.frame(matrix(NA, no.of.timesteps, 5))
colnames(sim.table) <- c('WeekNo','S','I1','I2','R')
tmp <- NA

increment <- function(pop,rno,lambda){
  if(pop==0){
    if(rno < lambda){
      pop <- pop+1
    }
  } else if(pop==1 | pop==2){
    pop <- pop+1
  }
  return(pop)
}

time_increment <- function(x){
  lambda <- beta*sum(x==1| x==2)/length(x)
  random.no <- runif(length(x))
  
  tmp <- c(sum(x==0),sum(x==1), sum(x==2), sum(x==3))
  pop <<- mapply(increment,pop=x, rno=random.no, lambda= lambda)
  
  tmp
}

for(i in 1:no.of.timesteps){
  sim.table[i,] <- c(i,(time_increment(pop)))
  
}

sim.table
