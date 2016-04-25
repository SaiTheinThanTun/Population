####comparing ODE and IBM####
setwd("C:/wd")
c.names <- c('no','time','Sh','X','lam_h','S','Z','lam') #setting column names

ode <- read.csv('summary_ode_.csv')
#ode <- read.csv('avg_summary_ibm_original2.csv')
#remove the last Y column
ode <- ode[,-9]
colnames(ode) <- c.names

ibm <- read.csv('avg_summary_ibm_2016-04-25.csv.csv')
colnames(ibm) <- c.names


#head(ode)
#head(ibm)


#difference of ode and ibm
dif <- ode-ibm
hist(dif[,3])
plot(dif[,3],dif[,3])
plot(1:nrow(dif),dif[,3])

#difference between lambda_h
plot(1:nrow(dif),dif[,5])
