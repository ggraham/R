annePlot<-function(melted_data){
  a<-ggplot(data = melted_data, aes(x = n_exons, y = n_trans, colour = L1))
#  a<-a+geom_point(pch = 19, na.rm = TRUE, alpha=0.2)
  a<-a+geom_smooth(method="loess")
  a<-a+xlim(1, 50)
  a<-a+ylim(1, 25)
#  a<-a+scale_y_log10()
#  a<-a+scale_x_log10()
#  a<-a+geom_density2d()
#  a<-a+geom_bar(stat="identity")
  return(a)
}