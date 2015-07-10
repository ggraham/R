latticePSI_RI<-function(identifier=NULL, ...){
  print(
    barchart(fac.order~as.vector(miso_RI_data[identifier,25:27]),
             xlim=c(0, 1),
             col=c(rgb(204,0,17,255,maxColorValue = 255),
                   rgb(51,204,255,255,maxColorValue = 255),
                   rgb(255,136,0,255,maxColorValue = 255)),
             stack=FALSE,
             box.ratio=20000,
             lwd=2,
             par.settings=list(fontsize=list(text=12,points=10)),
             xlab="Retained Intron PSI"),
    newpage=FALSE,
    prefix="plot")
}

latticePSI_SE<-function(identifier=NULL, ...){
  print(
    barchart(fac.order~as.vector(miso_SE_data[identifier,25:27]),
             xlim=c(0, 1),
             col=c(rgb(204,0,17,255,maxColorValue = 255),
                   rgb(51,204,255,255,maxColorValue = 255),
                   rgb(255,136,0,255,maxColorValue = 255)),
             stack=FALSE,
             box.ratio=20000,
             lwd=2,
             par.settings=list(fontsize=list(text=12,points=10)),
             xlab="Skipped Exon PSI"),
    newpage=FALSE,
    prefix="plot")
}