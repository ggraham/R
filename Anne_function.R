require(rtracklayer)
require(ggplot2)
require(reshape2)
require(stringr)
#
annePlot<-function(melted_data, x = n_exons, y = n_trans){
  a<-ggplot(data = melted_data, aes(x = n_exons, y = n_trans, colour = L1))
#  a<-a+geom_jitter()
#  a<-a+geom_point(pch = 19, na.rm = TRUE, alpha=0.5)
  a<-a+geom_smooth()
  a<-a+xlim(1, 100)
  a<-a+ylim(1, 200)
  a<-a+coord_cartesian(ylim=c(1, 15))
#  a<-a+scale_y_log10()
  a<-a+scale_x_log10()
#  a<-a+geom_density2d(h=c(100, 2))
#  a<-a+geom_bar(stat="identity")
  return(a)
}
#
trans_N_Exons_N_Trans<-function(GFF_file, id_index="default"){
  GR<-import.gff2(GFF_file)
  if (id_index=="default"){
    id_index<-GR$gene_id
  } else if (id_index=="PB") {
    id_index<-unlist(str_extract_all(GR$gene_id, "PB.\\d+"))
  }
  GR$gene_id
  t<-tapply(X=GR, INDEX = id_index, FUN = function(x) sum(x$type=="exon"))
  n<-tapply(X=GR, INDEX = id_index, FUN = function(x) length(unique(x$transcript_id)))
#  k<-length(gene_ids)
#  bar<-txtProgressBar(min=1, max=k, initial=1, style = 3)
#  for ( i in 1:k ){
#    m<-c(m, abs(min(start(GR[which(as.character(GR$gene_id) == gene_ids[i])])) - 
#                  max(end(GR[which(as.character(GR$gene_id) == gene_ids[i])]))))
#    t<-c(t, sum(GR[GR$gene_id==gene_ids[i]]$type=="exon"))
#    n<-c(n, sum(as.character(GR$gene_id) == gene_ids[i]))
#    setTxtProgressBar(bar, i)
#  }
#  close(bar)
  return(data.frame(n_exons=t, n_trans=n))
}