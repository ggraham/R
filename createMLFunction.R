require(caret)
require(doMC)
registerDoMC(8)

createMLFunction<-function(data_PSI, data_FPKMs, index){

  filter_vec_FPKM<-apply(data_FPKMs, MARGIN = 1, FUN = "var")
  FPKMs_filter<-as.matrix(data_FPKMs[(filter_vec_FPKM>=quantile(filter_vec_FPKM,0.25)),])
  FPKMs_filter<-FPKMs_filter[apply(FPKMs_filter, MARGIN = 1, FUN = function(x) median(x)>=0.5),]
  
  filter_vec_PSI<-apply(data_PSI, MARGIN = 1, FUN = "var")
  PSIs_filter<-as.matrix(data_PSI[(filter_vec_PSI>=quantile(filter_vec_PSI,0.1)),])
  PSIs_filter<-PSIs_filter[apply(PSIs_filter, MARGIN = 1, FUN = function(x) all(any(x>=0.1), any(x<=0.9))),]
  
  if(missing(index)){
    event<-1:nrow(PSIs_filter)
  } else {
    event<-index
  }
  print("Low variance variables removed")
  print(paste(nrow(FPKMs_filter), "variable rows"))
  print(paste(nrow(PSIs_filter), "response vector rows"))
  
  #grid<-expand.grid(n.trees=c(5000),
  #                  interaction.depth=c(1,2),
  #                  shrinkage=c(0.1,1),
  #                  n.minobsinnode=c(10))
  #grid<-expand.grid(alpha=c(0.001, 0.01, 0.05, 0.1), lambda=c(0.00001, 2, 10))
  grid<-expand.grid(C=c(0.00001), sigma=c(0.05))
  
  trC<-trainControl(method = "repeatedcv", 
                    number = 5,
                    repeats = 2,
                    #timingSamps=1,
                    verboseIter = FALSE, 
                    allowParallel = TRUE,
                    classProbs=FALSE)
  
  set.seed(5)
  model_list<-list()
  TIME<-Sys.time()
  for (i in event){
    cat("\r",i)
    resp<-PSIs_filter[i,]
    model_list[[i]]<-train(x = t(FPKMs_filter), 
                           y = as.numeric(resp),
                           method = "svmRadial",
                           preProcess = c("center", "scale"),
                           tuneGrid=grid,
                           trControl = trC,
                           metric="RMSE")
    
  }
  cat("\n", Sys.time() - TIME)
  return(model_list)
}