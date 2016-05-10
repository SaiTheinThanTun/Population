
###generic function for age distribution####
ex_age <- function(filename='age_union.csv'){
  age <- read.csv(filename, header=FALSE)[,1] #age vs population
  age <- read.csv('age_union.csv', header=FALSE)[,1] 
  #this can be automated to read in several file name as list and use lapply
  l.age <- length(age) #length of age vector##also the index for the last age group which would be 98+ or 90+, etc
  total.age <- sum(age)
  l2model <- round(l.age*.1) #length to model #start from 20% back from last
  
  data.age <- tail(age, l2model)
  data.age.target <- sum(data.age) #target population to extrapolate
  data.age <- data.age[-(length(data.age))]
  
  age_fun <- function(g){
    model.age <<- tmp <- NA
    model.age[1] <<- data.age[1]
    for(i in 1:(l2model-1)){
      tmp[i] <- data.age[i]-model.age[i]
      model.age[i+1] <<- model.age[i]*exp(g) #model.age is up to the last aggregated group. ie l2model
    }
    sum(tmp^2)
  }
  m <- optimise(age_fun, interval=c(-3,3))$minimum #value to extrapolate the age str
  
  
  
  ##while loop####
  j <- l2model
  while(sum(model.age)<data.age.target){
    model.age[j+1] <- model.age[j]*exp(m)
    j <- j+1
  }
  
  if(sum(model.age)!=data.age.target){ #being careful of the unexpected case
    diff <- sum(model.age)-data.age.target #put the reminder of the population to the last value of the data 
    if(model.age[length(model.age)]-diff >0){
      model.age[length(model.age)] <- model.age[length(model.age)]-diff
    } else 
    {
      model.age[length(model.age)+1] <- diff 
    }
  }
  newage <- append(age[1:(length(age)-l2model)], model.age) #, after= length(age)-l2model)
  if(sum(newage)==total.age){
    write.csv(newage, file=paste('ex_age_',filename, sep=''))
  } else { print("Simulated total doesn't match with the real total!")}
}

#compareAge <- cbind(age,newage)
##draft####
# age-exp <- function(x){
#   y <- m*x + b
# }
# 
# smallest <- function(x){
#   abs(sin(730*x)*.02)
# }
# 
# 
# optimize(smallest, interval=c(0.005,0.01))
# 
# #extrapolating the age structure over 80
# age_prob_0to97 <- read.csv("C:/wd/0to97_age_prob.csv", header=FALSE)
# age <- 0:97
# age_after80 <- tail(age, 18)
# y <- tail(age_prob_0to97[,1], 18)


# ####Fitting age structure with optimize####
# pop_after80 <- read.csv("C:/wd/age struct/actual_age_pop_aft80.csv", header=FALSE)[,1]
# pop_after80[1] <- 90000
# tmp <- NA
# age_fun <- function(g){
#   popMod <<- tmp <- NA
#   popMod[1] <<- 90000
#   for(i in 1:(length(pop_after80)-1)){
#     tmp[i] <- pop_after80[i]-popMod[i]
#     popMod[i+1] <<- popMod[i]*exp(g)
#     #tmp[i] <- pop_after80[i+1]-(pop_after80[i]*exp(g))
#   }
#   sum(tmp^2)
# }
# optimise(age_fun, interval=c(-3,3))
# 
# #the result is -0.1680289
# #Lisa's result from fitting with the eyes in Excel is -.19
# 
# age_fun(-0.1680289)
# a.17 <- popMod
# age_fun(-.19)
# b.19 <- popMod
# 
# 
# plot(pop_after80, type= 'p', col= 'black')
# lines(a.17, col= 'red')
# lines(b.19, col='blue')
# 
# compareModels <- cbind(pop_after80,a.17,b.19)
# 
# popMod <<- tmp <- NA
# popMod[1] <- 90000
# uptoage <- 101-80 #101
# for(i in 1:uptoage){
#   tmp[i] <- pop_after80[i]-popMod[i]
#   popMod[i+1] <- popMod[i]*exp(-0.1680289)
# }
# #98+ year is 7967
# sum(tail(popMod,4))

