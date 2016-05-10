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


####Fitting age structure with optimize####
pop_after80 <- read.csv("C:/wd/age struct/actual_age_pop_aft80.csv", header=FALSE)[,1]
pop_after80[1] <- 90000
tmp <- NA
age_fun <- function(g){
  popMod <<- tmp <- NA
  popMod[1] <<- 90000
  for(i in 1:(length(pop_after80)-1)){
    tmp[i] <- pop_after80[i]-popMod[i]
    popMod[i+1] <<- popMod[i]*exp(g)
    #tmp[i] <- pop_after80[i+1]-(pop_after80[i]*exp(g))
  }
  sum(tmp^2)
}
optimise(age_fun, interval=c(-3,3))

#the result is -0.1680289
#Lisa's result from fitting with the eyes in Excel is -.19

age_fun(-0.1680289)
a.17 <- popMod
age_fun(-.19)
b.19 <- popMod


plot(pop_after80, type= 'p', col= 'black')
lines(a.17, col= 'red')
lines(b.19, col='blue')

compareModels <- cbind(pop_after80,a.17,b.19)
