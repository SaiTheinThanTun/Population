
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
    write.table(newage, file=paste('ex_age_',filename, sep=''), row.names=FALSE, col.names=FALSE, sep=',')
  } else { print("Simulated total doesn't match with the real total!")}
}

