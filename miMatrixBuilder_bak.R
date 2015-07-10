miMatrixBuilder<-function(data_PSIs, data_FPKMs, indices=1:ncol(data_PSIs),i2=NULL, i1=NULL, boot=1, R=10){

  require(entropy)
  #filtering
  filter_vec_FPKM<-apply(data_FPKMs, MARGIN = 1, FUN = "var")
  FPKMs_filter<-as.matrix(data_FPKMs[(filter_vec_FPKM>=quantile(filter_vec_FPKM,0.25)),])
  FPKMs_filter<-FPKMs_filter[apply(FPKMs_filter, MARGIN = 1, FUN = function(x) median(x)>=0.5),]
  #print(dim(FPKMs_filter))
  filter_vec_PSI<-apply(data_PSIs, MARGIN = 1, FUN = "var")
  PSIs_filter<-as.matrix(data_PSIs[(filter_vec_PSI>=quantile(filter_vec_PSI,0.1)),])
  PSIs_filter<-PSIs_filter[apply(PSIs_filter, MARGIN = 1, FUN = function(x) all(any(x>=0.1), any(x<=0.9))),]
  #print(dim(PSIs_filter))

  
  if (is.null(i1)){
    i1<-1:nrow(PSIs_filter)
  }
  if (is.null(i2)){
    i2<-1:nrow(FPKMs_filter)
  }
  out_mat<-matrix(ncol = length(i1), nrow = length(i2))
  out_mat_error<-matrix(ncol = length(i1), nrow = length(i2))
  error_var<-numeric(length = R)
  #event loop
  for ( i in i1){
#    miso_disc<-discretize(as.numeric(PSIs_filter[i,]), numBins = 5)
  #  out_mat[,i]<-foreach( n = i2, .combine=c) %dopar% mi.Dirichlet(cbind(miso_disc, discretize(as.numeric(FPKMs_filter[n,]), numBins=5)), a=1/6)
    for ( n in i2){
 #     fpkm_disc<-discretize(as.numeric(FPKMs_filter[n,]), numBins = 5)
 #     out_mat[(n+1-(min(i2))),i]<-mi.Dirichlet(cbind(miso_disc, fpkm_disc), a = 1/6)
      for ( k in 1:R ){
        #bootstrapping
        indices<-sample(24, replace = FALSE, size = floor(ncol(PSIs_filter)*boot))
#        print(indices)
        miso_disc_boot<-discretize(as.numeric(PSIs_filter[i,indices]), numBins = 5)
        fpkm_disc_boot<-discretize(as.numeric(FPKMs_filter[n,indices]), numBins = 5)
        error_var[k]<-mi.Dirichlet(discretize2d(x1 = as.numeric(PSIs_filter[i,indices]), x2 = as.numeric(FPKMs_filter[n,indices]), numBins1 = 5, numBins2 = 5), a = 1/25)
      }
      #print(error_var)
      out_mat_error[(n+1-(min(i2))),i]<-sd(error_var)
      out_mat[(n+1-(min(i2))),i]<-mean(error_var)
      
    }
  }
  row.names(out_mat)<-row.names(FPKMs_filter)[i2]
  row.names(out_mat_error)<-row.names(FPKMs_filter)[i2]
  colnames(out_mat)<-row.names(PSIs_filter)[i1]
  colnames(out_mat_error)<-row.names(PSIs_filter)[i1]
  return(list(MI=out_mat, MI_SD=out_mat_error))
  
}