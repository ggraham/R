miListInput<-function(data_PSIs, data_FPKMs, R){
  filter_vec_FPKM<-apply(data_FPKMs, MARGIN = 1, FUN = "var")
  FPKMs_filter<-as.matrix(data_FPKMs[(filter_vec_FPKM>=quantile(filter_vec_FPKM,0.25)),])
  FPKMs_filter<-FPKMs_filter[apply(FPKMs_filter, MARGIN = 1, FUN = function(x) median(x)>=0.5),]
  print(dim(FPKMs_filter))
  filter_breaks_FPKM<-apply(FPKMs_filter, MARGIN=1, FUN = function(x) length(unique(x))>7)
  FPKMs_filter<-FPKMs_filter[filter_breaks_FPKM,]
  FPKMs_filter<-as.matrix(FPKMs_filter)
  print(dim(FPKMs_filter))
  #print(dim(FPKMs_filter))
  filter_vec_PSI<-apply(data_PSIs, MARGIN = 1, FUN = "var")
  PSIs_filter<-as.matrix(data_PSIs[(filter_vec_PSI>=quantile(filter_vec_PSI,0.25)),])
  PSIs_filter<-PSIs_filter[apply(PSIs_filter, MARGIN = 1, FUN = function(x) all(any(x>=0.1), any(x<=0.9))),]
  print(dim(PSIs_filter))
  filter_breaks_PSI<-apply(PSIs_filter, MARGIN=1, FUN = function(x) length(unique(x))>7)
  PSIs_filter<-PSIs_filter[filter_breaks_PSI,]
  PSIs_filter<-as.matrix(PSIs_filter)
  print(dim(PSIs_filter))
  #
  #gr<-expand.grid(psi.n=row.names(PSIs_filter), gene.n=row.names(FPKMs_filter))
  #gr.dat<-cbind(PSIs_filter[gr$psi.n,], FPKMs_filter[gr$gene.n,])
  #
  #list_out<-vector("list", length = nrow(PSIs_filter))
  #for (i in 1:nrow(PSIs_filter)){
  #  list_out[[i]]<-vector("list", length = nrow(FPKMs_filter))
  #  for (n in 1:nrow(FPKMs_filter)){
  #    list_out[[i]][[n]]<-vector("list", length = R)
  #    for (k in 1:R){
  #      ind<-sample(24, 24, replace = TRUE)
  #      list_out[[i]][[n]][[k]]<-discretize2d(x1 = PSIs_filter[i,ind], x2 = FPKMs_filter[n,ind], numBins1 = 5, numBins2 = 5)
  #    }
  #  }
  #}
  return()
}