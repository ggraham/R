clusterRunFunc<-function(cl){
  require(snow)
  matrix_out_list<-clusterApply(cl,x = list(i2=1:2039, 
                                            i2=2040:4079,
                                            i2=4080:6120,
                                            i2=6121:8159, 
                                            i2=8160:10199,
                                            i2=10200:12236), fun = miMatrixBuilder, data_PSIs=miso_RI_data[,1:24], data_FPKMs=gene_FPKMs[,1:24], boot=0.75, R=100, i1=1:5)
  stopCluster(cl)
  return(matrix_out_list)
}
combineFunction<-function(listed){
  out_mat<-listed[[1]]
  for ( i in 2:(length(listed)-1)){
    out_mat<-rbind(out_mat, listed[[i]])
  }
  return(out_mat)
}