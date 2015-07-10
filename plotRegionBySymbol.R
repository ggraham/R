require(biomaRt)
require(Gviz)
require(lattice)
makeTracks<-function(events=NULL, dat=dat){
  if(length(events)!=0){
    RI_AT<-AnnotationTrack(range = miso_RI_data_AnnTrack_GR[events[[2]]],
                           from = dat$start_position, 
                           to = dat$end_position, 
                           chromosome = paste("chr", dat$chromosome_name, sep=""),
                           genome = "hg19", 
                           id = row.names(miso_RI_data[events[[2]],]), 
                           name="Retained Intron", 
                           fun=latticePSI_RI,
                           details.size=0.85)
    SE_AT<-AnnotationTrack(range = miso_SE_data_AnnData_GR[events[[1]]],
                           from = dat$start_position, 
                           to = dat$end_position, 
                           chromosome = paste("chr", dat$chromosome_name, sep=""),
                           genome = "hg19", 
                           id = row.names(miso_SE_data[events[[1]],]), 
                           name="Skipped Exon", 
                           fun=latticePSI_SE,
                           details.size=0.85)
  } else {
    RI_AT<-AnnotationTrack(range = miso_RI_data_AnnTrack_GR,
                           from = dat$start_position, 
                           to = dat$end_position, 
                           chromosome = paste("chr", dat$chromosome_name, sep=""),
                           genome = "hg19", 
                           id = row.names(miso_RI_data), 
                           name="Retained Intron", 
                           fun=latticePSI_RI,
                           details.size=0.85)
    SE_AT<-AnnotationTrack(range = miso_SE_data_AnnData_GR,
                           from = dat$start_position, 
                           to = dat$end_position, 
                           chromosome = paste("chr", dat$chromosome_name, sep=""),
                           genome = "hg19", 
                           id = row.names(miso_SE_data), 
                           name="Skipped Exon", 
                           fun=latticePSI_SE,
                           details.size=0.85)
  }
  return(list(SE_AT=SE_AT, RI_AT=RI_AT))
}
checkEvents<-function(coord_dat, track_list){
  track_out<-c(hg19_GR)
  coord_dat_GRanges<-GRanges(seqnames = paste("chr", coord_dat$chromosome_name, sep=""),
                             ranges = IRanges(start = coord_dat$start_position,
                                              end = coord_dat$end_position,
                                              names = coord_dat$hgnc_symbol),
                             strand = "*")
  if(sum(countOverlaps(coord_dat_GRanges, miso_SE_data_AnnData_GR))!=0){
    track_out<-c(track_out, track_list[[2]])
  }
  if(sum(countOverlaps(coord_dat_GRanges, miso_RI_data_AnnTrack_GR))!=0){
    track_out<-c(track_out, track_list[[3]])
  }
  return(c(track_out, track_list[[4]]))
}
  
plotRegionBySymbol<-function(sy=NULL,events=NULL){
  biom<-useDataset(dataset = "hsapiens_gene_ensembl",
                   mart = useMart(host="feb2014.archive.ensembl.org",
                                  biomart="ENSEMBL_MART_ENSEMBL"))
  dat<-getBM(attributes = c("start_position", "end_position", "chromosome_name", "hgnc_symbol"),
        filters = c("hgnc_symbol"), mart=biom, values=sy)
  dat$start_position<-dat$start_position-1500
  dat$end_position<-dat$end_position+1500
  print(dat)
  dat<-dat[nchar(as.character(dat$chromosome_name))==min(nchar(as.character(dat$chromosome_name))),]
  print(dat)
  eventAnnotation<-makeTracks(events, dat)
  #eventAnnotation$SE_AT
  #print(eventAnnotation$SE_AT)
  trackvec<-c(hg19_GR, eventAnnotation$SE_AT, eventAnnotation$RI_AT, tx.hq)
  print(length(trackvec))
  trackvec<-checkEvents(dat, trackvec)
  print(length(trackvec))
  sizes<-c(1, 2)
  if (length(trackvec)==3){
    sizes<-c(1, 3, 2)
  } else if (length(trackvec)==4){
    sizes<-c(1, 3, 3, 2)
  }
  #trackvec<-c(hg19_GR, tx.hq)
  for (i in 1:length(trackvec)){
    displayPars(trackvec[[i]])<-list(fontcolor.title="black", cex.title=0.9, fill="black")
  }
  #plotTracks(trackvec[c(1,3,4)], 
  #           sizes=c(1,3,3),
  plotTracks(trackvec,
             sizes=sizes,
             from = dat$start_position, 
             to = dat$end_position, 
             chromosome = paste("chr", dat$chromosome_name, sep=""), 
             main = dat$hgnc_symbol, 
             transcriptAnnotation="symbol"
             #type=c("p", "histogram", "g"), 
             #showSampleNames=TRUE, 
             #groups=c("TC32 shFli1", "TC32 YK"), 
             #ylim=c(-1, 1)
  )
}