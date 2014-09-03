
pwmToCoverage <- function(pwm,genome,min="70%",removeRand=FALSE){
  
  allchrs <- seqnames(genome)
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  intergerList <- lapply(allchrs,function(x)pwmHitAsCoverage(pwm,genome,min,x))
  myrle <- RleList(intergerList,compress=F)
  names(myrle) <- allchrs
  myrle
}  

pwmHitAsCoverage <- function(pwm,genome,min,chrofinterest){
  posMotifs <- matchPWM(pwm,genome[[chrofinterest]],min.score=min)
  negMotifs <- matchPWM(reverseComplement(pwm),genome[[chrofinterest]],min.score=min)
  if(length(posMotifs) > 0){
    rleMotifHitPos <- coverage(GRanges(seqnames=chrofinterest,ranges(posMotifs),strand="+"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitPos <- RleList(rep(0,length(genome[[chrofinterest]])))
  }
  if(length(negMotifs) > 0){
    rleMotifHitNeg <- coverage(GRanges(seqnames=chrofinterest,ranges(negMotifs),strand="-"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitNeg <- RleList(rep(0,length(genome[[chrofinterest]])))    
  }
  rleTotal <- rleMotifHitPos+rleMotifHitNeg
  return(rleTotal[[1]])
}


motifCov <- function(genome,regions,pwm,chrOfInterest){
  regionViews <- Views(genome[[chrOfInterest]],ranges(regions[seqnames(regions) %in% chrOfInterest]))
  trial <- matchPWM(pwm,regionViews,min.score = 0,with.score = T)
  #theRanges <- resize(as(trial,"IRanges"),1,"center")
  theRanges <- as(trial,"IRanges")
  if(length(theRanges) > 0){
    motifCov <- unlist(coverage(GRanges(chrOfInterest,theRanges,"*",elementMetadata(trial)$score),weight="elementMetadata.trial..score"))
  }else{
    motifCov <- unlist(RleList(rep(0,length(genome[[chrOfInterest]]))))
  }
  return(motifCov)
}

makeMotifScoreRle <- function(pwm,regions,genome,extend,removeRand=FALSE,strandScore="mean"){
  regions <- GRanges(seqnames(regions),IRanges(start(regions)-extend,end(regions)+extend),strand=strand(regions),elementMetadata(regions))
  lengths <- seqlengths(genome)
  ## Filter testRanges to those contained within chromosomes.
  message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
  testRangeNames <- unique(seqnames(regions))
  temptestranges <- GRanges()
  maxDistance <- extend
  for(i in 1:length(testRangeNames)){
    perchrRanges <- regions[seqnames(regions) %in% as.vector(testRangeNames[i])]
    temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                    & start(perchrRanges)-maxDistance > 0 ])
    #print(i)
  }
  
  message("..Done")
  message("Filtered ",length(regions)-length(temptestranges)," of ",length(regions)," regions")
  
  regions <- temptestranges
  
  allchrs <- as.vector(unique(seqnames(regions)))
  
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  forMotif <- list()
  revMotif <- list()
  message("Scoring motifs on positive strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])
    forMotif[[k]] <- unlist(motifCov(genome,regions,pwm,allchrs[k]))  
  }
  message("..done")
  
  message("Scoring motifs on negative strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])    
    revMotif[[k]] <- unlist(motifCov(genome,regions,reverseComplement(pwm),allchrs[k]))    
  }
  message("..done")
  
  #motifScoreRLEComplement <- lapply(allchrs,function(x)motifCov(genome,regions,reverseComplement(pwm),x))
  #myrle <- RleList(motifScoreRLE,compress=F)
  # myrleComplement <- RleList(motifScoreRLEComplement,compress=F)
  revMotif <- RleList(revMotif,compress=F)
  forMotif <- RleList(forMotif,compress=F)
  if(strandScore=="sum"){
    MotifScore <- revMotif+forMotif
  }
  if(strandScore=="mean"){  
    MotifScore <- (revMotif+forMotif)/2
  }
  if(strandScore=="max"){  
    MotifScore <- pmax(revMotif,forMotif)
  }  
  names(MotifScore) <- allchrs
  return(MotifScore)
}


library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
library(MotifDb)
library(Biostrings)
library(seqLogo)
#load("/home/pgellert/MatthiasTrial/superenhancers.RData")
load("/Users/tcarroll/Downloads/superenhancers.RData")

dpctcfPeaks <- ChIPQC:::GetGRanges("/home//pgellert/Dropbox (Lymphocyte_Developme)/tracktables/DP_CTCF_WithInput_DP_Input_peaks.bed")
#dpctcfPeaks <- ChIPQC:::GetGRanges("/Users/tcarroll/Downloads/DP_CTCF_WithInput_DP_Input_peaks.bed")


mdb.ctcf <- MotifDb [grep ('ctcf', values (MotifDb)$geneSymbol, ignore.case=TRUE)]
#mdb.ctcf[[2]]

png("/home/pgellert/MatthiasTrial/seqLogoMotif1.png")
seqLogo(mdb.ctcf[[1]])
dev.off()

png("/home/pgellert/MatthiasTrial/seqLogoMotif2.png")
seqLogo(mdb.ctcf[[2]])
dev.off()


motif70 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="70%",removeRand=T)
motif75 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="75%",removeRand=T)
motif80 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="80%",removeRand=T)
motif85 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="85%",removeRand=T)
motif90 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="90%",removeRand=T)


ctcfDP <- regionPlot("/home/pgellert/MatthiasTrial/DP_thymocyte_CTCF_Shih_sorted.bam",DP_thymocytes,nOfWindows=1000,style="region",format="bam",FragmentLength=130)
cohesinDP <- regionPlot("/home/pgellert/MatthiasTrial/DP_thymocyte_Rad21_all_m1_v2_Vlad_sorted.bam",DP_thymocytes,nOfWindows=1000,style="region",format="bam",FragmentLength=130)


ctcfMotif80 <- regionPlot(motif80,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)

ctcfMotifMeanScoreBnary80 <- regionPlot(motif80,dpctcfPeaks,style="point",format="rlelist",FragmentLength=130)

ctcfMotif85 <- regionPlot(motif85,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)

ctcfMotifMeanScoreBnary85 <- regionPlot(motif80,dpctcfPeaks,style="point",format="rlelist",FragmentLength=130)

gc()

rleMotifScores1 <- makeMotifScoreRle(mdb.ctcf[[1]],dpctcfPeaks,Mmusculus,1500,removeRand=TRUE,strandScore="mean")
gc()
ctcfMotifMeanScore <- regionPlot(rleMotifScores1,dpctcfPeaks,style="point",format="rlelist",FragmentLength=130)




genome <- Mmusculus
chrOfInterest <- "chr19"
dpctcfPeaks <- ChIPQC:::GetGRanges("/home//pgellert/Dropbox (Lymphocyte_Developme)/tracktables/DP_CTCF_WithInput_DP_Input_peaks.bed")
chr19Peaks <- dpctcfPeaks[seqnames(dpctcfPeaks) %in% "chr19"]
chr19Peaks <- chr19Peaks[-length(chr19Peaks)]
rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],chr19Peaks,2000,removeRand=FALSE,strandScore="mean")
ctcfMotifMeanScore <- regionPlot(rleMotif1,chr19Peaks,style="point",format="rlelist",FragmentLength=130)

chr1Peaks <- dpctcfPeaks[seqnames(dpctcfPeaks) %in% "chr1"]
chr1Peaks <- chr1Peaks[-length(chr1Peaks)]
rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],chr1Peaks,2000,removeRand=FALSE,strandScore="mean")
ctcfMotifMeanScore <- regionPlot(rleMotif1,chr1Peaks,style="point",format="rlelist",FragmentLength=130)

Top5000 <- dpctcfPeaks[order(as.numeric(as.vector(dpctcfPeaks$Score)),decreasing=T)][1:1000]
rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],Top5000,Mmusculus,2000,removeRand=TRUE,strandScore="mean")
ctcfMotifMeanScore <- regionPlot(rleMotif1,Top5000,style="point",format="rlelist",FragmentLength=130)

names(MotifSumScore) <- allchrs
names(MotifMeanScore) <- allchrs

ctcfMotifSumScore <- regionPlot(MotifSumScore,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)
ctcfMotifMeanScore <- regionPlot(MotifMeanScore,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)


motifScoreRLE2 <- lapply(allchrs,function(x)motifCov(genome,regions,mdb.ctcf[[2]],x))
motifScoreRLEComplement2 <- lapply(allchrs,function(x)motifCov(genome,regions,reverseComplement(mdb.ctcf[[2]]),x))
myrle2 <- RleList(motifScoreRLE2,compress=F)
myrleComplement2 <- RleList(motifScoreRLEComplement2,compress=F)

MotifSumScore2 <- myrle2+myrleComplement2
MotifMeanScore2 <- (myrle2+myrleComplement2)/2


names(MotifSumScore2) <- allchrs
names(MotifMeanScore2) <- allchrs

ctcfMotifSumScore2 <- regionPlot(MotifSumScore2,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)
ctcfMotifMeanScore2 <- regionPlot(MotifMeanScore2,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)


rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],DP_thymocytes,10000,removeRand=FALSE,strandScore="mean")
ctcfMotifMeanScore <- regionPlot(rleMotif1,DP_thymocytes,nOfWindows=1000,style="region",format="rlelist",FragmentLength=130)


load("/home//pgellert/MatthiasTrial/Motif85.RData")

ctcfMotifMeanScoreBnary <- regionPlot(motif85,chr1Peaks,style="point",format="rlelist",FragmentLength=130)

plot(c(seq(1,3001),3001+seq(1,10000,100),10301+seq(1,3001)),runmean(10^colMeans(assays(ctcfMotifScore)[[1]]),50),type="l")


Top5000 <- dpctcfPeaks[order(as.numeric(as.vector(dpctcfPeaks$Score)),decreasing=T)][1:1000]
rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],dpctcfPeaks,Mmusculus,2000,removeRand=TRUE,strandScore="mean")
ctcfMotifMeanScore <- regionPlot(rleMotif1,dpctcfPeaks,style="point",format="rlelist",FragmentLength=130)
exptData(ctcfMotifMeanScore) <- list(names=c("sampleName"))
plotRegion(ctcfMotifMeanScore)

rleMotif1 <- makeMotifScoreRle(mdb.ctcf[[1]],DP_thymocytes,Mmusculus,2000,removeRand=TRUE,strandScore="mean")
ctcfMotifEnhMeanScore <- regionPlot(rleMotif1,DP_thymocytes,style="region",format="rlelist",FragmentLength=130)
exptData(ctcfMotifEnhMeanScore) <- list(names=c("sampleName"))
plotRegion(ctcfMotifEnhMeanScore)

MoreSets <- c(ChIPQC:::GetGRanges(DP_thymocytes,simplify=T),ChIPQC:::GetGRanges(ESCs,simplify=T))

rleMotif12 <- makeMotifScoreRle(mdb.ctcf[[1]],MoreSets,Mmusculus,2000,removeRand=TRUE,strandScore="mean")
ctcfMotifEnhMeanScore2 <- regionPlot(rleMotif12,MoreSets,style="region",nOfWindows=1000,format="rlelist",FragmentLength=130)
exptData(ctcfMotifEnhMeanScore2) <- list(names=c("sampleName"))
plotRegion(ctcfMotifEnhMeanScore2)

rleMotif14 <- makeMotifScoreRle(mdb.ctcf[[1]],DP_thymocytes,Mmusculus,6000,removeRand=TRUE,strandScore="mean")

ctcfMotifEnhMeanScore4 <- regionPlot(rleMotif14,DP_thymocytes,style="region",format="rlelist",FragmentLength=130,
                                    distanceInRegionStart = 5000,
                                    distanceOutRegionStart = 5000, distanceInRegionEnd = 5000,
                                    distanceOutRegionEnd = 5000, 
                                    )
exptData(ctcfMotifEnhMeanScore4) <- list(names=c("sampleName"))

plot(caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]]),50),type="l")

> mean(width(DP_thymocytes)/100)
[1] 376.4326
> median(width(DP_thymocytes)/100)
[1] 309.19
> plot(j$data$xIndex,caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]]),300),type="l")
> plot(j$data$xIndex,caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]]),3),type="l")
> 
#plotRegion(ctcfMotifEnhMeanScore4)

extendedDP <-  GRanges(seqnames(DP_thymocytes),IRanges(start(DP_thymocytes)-(width(DP_thymocytes)+1000),end(DP_thymocytes)+(width(DP_thymocytes)+1000)),strand="+",elementMetadata(DP_thymocytes))
  
rleMotif14 <- makeMotifScoreRle(mdb.ctcf[[1]],extendedDP,Mmusculus,1000,removeRand=TRUE,strandScore="mean")

ctcfMotifEnhMeanScore4 <- regionPlot(rleMotif14,DP_thymocytes,style="region",format="rlelist",FragmentLength=130,distanceInRegionStart = 5000,
                                     distanceOutRegionStart = 5000, distanceInRegionEnd = 5000,
                                     distanceOutRegionEnd = 5000)

exptData(ctcfMotifEnhMeanScore4) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore4)
plot(plotRegionRes$data$xIndex,caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]]),2000),type="l")




extendedDP <-  GRanges(seqnames(DP_thymocytes),IRanges(start(DP_thymocytes)-(width(DP_thymocytes)+1000),end(DP_thymocytes)+(width(DP_thymocytes)+1000)),strand="+",elementMetadata(DP_thymocytes))
rleMotif14 <- makeMotifScoreRle(mdb.ctcf[[1]],extendedDP,Mmusculus,1000,removeRand=TRUE,strandScore="mean")
motifScores_DP_thy_Enh <-rleMotif14
save(motifScores_DP_thy_Enh,file="/home//pgellert/MatthiasTrial/MotifScores.RData") 
  
ctcfMotifEnhMeanScore4 <- regionPlot(rleMotif14,DP_thymocytes,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 2000,
                                     distanceOutRegionStart = 5000, distanceInRegionEnd = 2000,
                                     distanceOutRegionEnd = 5000)
exptData(ctcfMotifEnhMeanScore4) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore4)

plot(plotRegionRes$data$xIndex,caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]]),500),type="l")

firstPart <- caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]])[1:7001],300)
mid <- caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]])[7002:8001],50)
SecondPart <- caTools:::runmean(colMeans(assays(ctcfMotifEnhMeanScore4)[[1]])[8002:15002],300)

plot(plotRegionRes$data$xIndex,c(firstPart,mid,SecondPart),type="l",xaxt='n')
axis(side=1,at=plotRegionRes$scales$scales[[1]]$breaks,labels=plotRegionRes$scales$scales[[1]]$labels,las=2)




###########


extendedDP <-  GRanges(seqnames(DP_thymocytes),IRanges(start(DP_thymocytes)-(2*width(DP_thymocytes)),end(DP_thymocytes)+(2*width(DP_thymocytes))),strand="+",elementMetadata(DP_thymocytes))
rleMotif2 <- makeMotifScoreRle(mdb.ctcf[[1]],extendedDP,Mmusculus,1000,removeRand=TRUE,strandScore="max")
motifScores_DP_thy_Enh_Max <-rleMotif2
save(motifScores_DP_thy_Enh_Max,file="/home//pgellert/MatthiasTrial/MaxStrand_MotifScores.RData") 

ctcfMotifEnhMeanScore42 <- regionPlot(rleMotif2,DP_thymocytes,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 2000,
                                     distanceOutRegionStart = 10000, distanceInRegionEnd = 2000,
                                     distanceOutRegionEnd = 10000)
exptData(ctcfMotifEnhMeanScore42) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore42)
plotRegionRes

ctcfMotifEnhMeanScore42 <- regionPlot(rleMotif2,DP_thymocytes,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 5000,
                                      distanceOutRegionStart = 10000, distanceInRegionEnd = 5000,
                                      distanceOutRegionEnd = 10000)
exptData(ctcfMotifEnhMeanScore42) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore42)
plotRegionRes


names(revMotif) <- allchrs
names(forMotif) <- allchrs

ctcfMotifEnhMeanScoreRev <- regionPlot(revMotif,DP_thymocytes,style="region",format="rlelist",nOfWindows=100,FragmentLength=130)

ctcfMotifEnhMeanScoreFor <- regionPlot(forMotif,DP_thymocytes,style="region",format="rlelist",nOfWindows=100,FragmentLength=130)

exptData(ctcfMotifEnhMeanScoreRev) <- list(names=c("sampleName"))
exptData(ctcfMotifEnhMeanScoreFor) <- list(names=c("sampleName"))
temp <- pmax(revMotif,forMotif)

ctcfMotifEnhMeanScoreMax <- regionPlot(temp,DP_thymocytes,style="region",format="rlelist",nOfWindows=100,FragmentLength=130)
exptData(ctcfMotifEnhMeanScoreMax) <- list(names=c("sampleName"))
plotRegion(ctcfMotifEnhMeanScoreMax)

extendedDP <-  GRanges(seqnames(DP_thymocytes),IRanges(start(DP_thymocytes)-(2*width(DP_thymocytes)),end(DP_thymocytes)+(2*width(DP_thymocytes))),strand="+",elementMetadata(DP_thymocytes))
rleMotif28 <- makeMotifScoreRle(mdb.ctcf[[1]],extendedDP,Mmusculus,1000,removeRand=TRUE,strandScore="max")

ctcfMotifEnhMeanScore42 <- regionPlot(rleMotif28,DP_thymocytes,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 5000,
                                      distanceOutRegionStart = 10000, distanceInRegionEnd = 5000,
                                      distanceOutRegionEnd = 10000)

exptData(ctcfMotifEnhMeanScore42) <- list(names=c("sampleName"))
plotRegion(ctcfMotifEnhMeanScore42)
########


extendedDP <-  GRanges(seqnames(DP_thymocytes),IRanges(start(DP_thymocytes)-(2*width(DP_thymocytes)),end(DP_thymocytes)+(2*width(DP_thymocytes))),strand="+",elementMetadata(DP_thymocytes))
rleMotif2 <- makeMotifScoreRle(mdb.ctcf[[1]],extendedDP,Mmusculus,1000,removeRand=TRUE,strandScore="max")
motifScores_DP_thy_Enh_Max <-rleMotif2

ctcfMotifEnhMeanScore42 <- regionPlot(rleMotif2,DP_thymocytes,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 2000,
                                      distanceOutRegionStart = 10000, distanceInRegionEnd = 2000,
                                      distanceOutRegionEnd = 10000)
exptData(ctcfMotifEnhMeanScore42) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore42)
plotRegionRes




MoreSets <- c(ChIPQC:::GetGRanges(DP_thymocytes,simplify=T),ChIPQC:::GetGRanges(ESCs,simplify=T))
MoreSets2 <-  GRanges(seqnames(MoreSets),IRanges(start(MoreSets)-(2*width(MoreSets)),end(MoreSets)+(2*width(MoreSets))),strand="+",elementMetadata(MoreSets))

rleMotif2 <- makeMotifScoreRle(mdb.ctcf[[1]],MoreSets2,Mmusculus,1000,removeRand=TRUE,strandScore="max")
motifScores_DP_thy_Enh_Max <-rleMotif2

ctcfMotifEnhMeanScore42 <- regionPlot(rleMotif2,MoreSets,style="region",format="rlelist",nOfWindows=1000,FragmentLength=130,distanceInRegionStart = 2000,
                                      distanceOutRegionStart = 10000, distanceInRegionEnd = 2000,
                                      distanceOutRegionEnd = 10000)
exptData(ctcfMotifEnhMeanScore42) <- list(names=c("sampleName"))
plotRegionRes <- plotRegion(ctcfMotifEnhMeanScore42)
plotRegionRes

plot(plotRegionRes$data$xIndex,apply(assays(ctcfMotifEnhMeanScore42)[[1]],2,function(x)mean(x,trim=0.2)),type="l",xaxt="n")
axis(side=1,at=plotRegionRes$scales$scales[[1]]$breaks,labels=plotRegionRes$scales$scales[[1]]$labels,las=2)

plot(plotRegionRes$data$xIndex,apply(assays(ctcfMotifEnhMeanScore42)[[1]],2,function(x)median(x)),type="l",xaxt="n")
axis(side=1,at=plotRegionRes$scales$scales[[1]]$breaks,labels=plotRegionRes$scales$scales[[1]]$labels,las=2)

plot(plotRegionRes$data$xIndex,caTools::runmean(apply(assays(ctcfMotifEnhMeanScore42)[[1]],2,function(x)median(x)),200),type="l",xaxt="n")
axis(side=1,at=plotRegionRes$scales$scales[[1]]$breaks,labels=plotRegionRes$scales$scales[[1]]$labels,las=2)



motif80DPEnh <- regionPlot(motif80,DP_thymocytes,style="percentOfRegion",distanceAround=100,format="rlelist")
