parApplyMI<-function(x, R=100){
  mi<-replicate(R, singleParApplyMI(x))
  return(c(MI=mean(mi), SD=sd(mi)))
}
singleParApplyMI<-function(x){
  require(entropy)
  ind<-sample(24, 24, replace = TRUE)
  d<-discretize2d(x1 = x[ind], x2 = x[(ind+24)], numBins1 = 7, numBins2 = 7)
  return(mi.Dirichlet(d, a = 1/49))
}