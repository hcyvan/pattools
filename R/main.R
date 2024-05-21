source('./R/base.R')

library(GenomicRanges)
library(ComplexHeatmap)
##########################################################################################
dataRatio<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold3.bed"),header = TRUE)
dataDepth<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold3.depth.bed"),header = TRUE)
dataRatio<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold3.bed"),header = TRUE)
dataDepth<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold3.depth.bed"),header = TRUE)
dataRatio<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold1.bed"),header = TRUE)
dataDepth<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold1.depth.bed"),header = TRUE)
####################################################################################
ratio<-mmRatio[mmRatio!=-1]
saveImage("hist.base.depth1.methylation.pdf",width=4, height=3.5)
hist(ratio,breaks = 200,main="CpGs of All samples",xlab="Methylation (%)")
dev.off()

###################################################################################
clincal<-Clincal()
dataRatio<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold10.bed"),header = TRUE)
mm<-dataRatio[sample(1:nrow(dataRatio),1000),-(1:3)]
mm<-clincal$pickSubjectIdByDataID(mm)
#mm<-clincal$pickCellTypeByDataID(mm)
Heatmap(mm)

##################################################################################
dataRatio<-loadData(file.path(CONFIG$DataInter,"CAC.depthMinThreshold1.bed"),header = TRUE)
dataRatio[dataRatio==-1]<-NA
mm<-clincal$pickSubjectIdByDataID(dataRatio)
data<-cbind(dataRatio[,1:3],mm)
colnames(data)[1]<-'chrom'
gr2<-GRanges(data)

marker<-loadData(file.path(CONFIG$DataRaw,"Cancer12MarkerV4.csv"))
marker<-dplyr::select(marker, chrom,start, end, hyper, Gene, abbr)
gr1<-GRanges(marker)

overlaps<-findOverlaps(gr1, gr2)
overlapDf<-data.frame(overlaps)

aa<-sapply(split(overlapDf, overlapDf$queryHits), function(x){
  target<-gr2[x$subjectHits]
  colMeans(data.frame(target@elementMetadata@listData),na.rm = TRUE)
})%>%t

marker[as.numeric(rownames(aa)),]

