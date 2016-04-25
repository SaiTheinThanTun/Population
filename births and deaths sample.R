#1-(1-(a*b*m))^z
#1-(1-(a*b))^(z*m)

S <- 600
I <- 200
muo <- .05
mui <- .05

timesteps <- 56
out <-matrix(NA,timesteps, 2)
out <- rbind(c(S,I),out)
for(i in 1:timesteps+1){
  #deaths and births
  out[i,1] <- out[i-1,1]- (out[i-1,1]*muo) + (out[i-1,1]*mui)
  out[i,2] <- out[i-1,2]- (out[i-1,2]*muo) + (out[i-1,2]*mui)
}
