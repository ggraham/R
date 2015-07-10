require(caret)
require(doMC)
registerDoMC(10)

createIsoFunction<-function(data_PSI, data_FPKMs, index=1:nrow(data_PSI)){
  event<-c(index)
  filter_vec_FPKM<-apply(data_FPKMs, MARGIN = 1, FUN = "var")
  FPKMs_filter<-as.matrix(data_FPKMs[(filter_vec_FPKM>=quantile(filter_vec_FPKM,0.25)),])
  filter_vec_PSI<-apply(data_PSI, MARGIN = 1, FUN = "var")
  PSIs_filter<-as.matrix(data_PSI[(filter_vec_PSI>=quantile(filter_vec_PSI,0.1)),])
  
  print("Low variance variables removed")
  print(paste(nrow(FPKMs_filter), "variable rows"))
  print(paste(nrow(PSIs_filter), "response vector rows"))
  
  #grid<-expand.grid(n.trees=c(5000),
  #                  interaction.depth=c(1,2),
  #                  shrinkage=c(0.1,1),
  #                  n.minobsinnode=c(10))
  grid<-expand.grid(alpha=c(0.001, 0.01, 0.05, 0.1), lambda=c(0.00001, 2, 10))
  
  trC<-trainControl(method = "cv", 
                    number = 5,
                    repeats = 5,
                    #timingSamps=1,
                    verboseIter = FALSE, 
                    allowParallel = TRUE,
                    classProbs=FALSE)
  
  set.seed(5)
  model_list<-list()
  
  for (i in event){
    print(i)
    resp<-PSIs_filter[i,]
    model_list[[i]]<-train(x = t(FPKMs_filter), 
                           y = as.numeric(resp),
                           method = "glmnet",
                           preProcess = c("center", "scale"),
                           tuneGrid=grid,
                           trControl = trC,
                           metric="RMSE")
    
  }
  return(model_list)
}