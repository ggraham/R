require(ggplot2)
require(reshape2)

plotIsoformFPKMs<-function(PB_ID=NULL, data=isoform_FPKMs, pheno_factor=as.factor(colnames(isoform_FPKMs)[25:27])){
  subset_FPKMs<-data[grep("PB_ID", row.names(data)),25:27]
  for (i in 1:ncol(subset_FPKMs)){
    subset_list<-list(colnames(subset_FPKMs)[i]=subset_FPKMs[,i]
  }
  melted<-melt(subset_list)
  ggplot(melted, 
  
}