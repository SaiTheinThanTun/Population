## MAKE SQUARE GRID WITH INDICES ##
yy <- seq(12.661,13.138,by=0.0458)
xx <- seq(102.492,102.787,by=0.0472)
ind <- matrix(1:60,ncol=6)



#////////////Catch people move into each grid/////////////////

for(j in 1:nrow(m)){
  pointidxs <- ind[findInterval(y coordinate for each person, yy), 
                   findInterval(x coordinate for each person, xx)]
  patch<-pointidxs
}


xx <- 1:4
yy <- 1:4
ind <- matrix(1:16,ncol=4)

for(j in 1:nrow(m)){
  pointidxs <- ind[findInterval(3, yy), 
                   findInterval(2, xx)]
  patch<-pointidxs
}