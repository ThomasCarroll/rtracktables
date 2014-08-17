
findconsensusRegions <- function(testRanges,bamFiles=NULL,method="majority",summit="mean",resizepeak="asw",overlap="any",fragmentLength=NULL,
                                 NonPrimaryPeaks=list(withinsample="drop",betweensample="mean")){

  testRanges <- GRangesList(
    bplapply(
      testRanges,
      ChIPQC:::GetGRanges)
    )
  
  ans <- lapply(
    names(bamFiles),
    function(x) summitPipeline(
      unlist(bamFiles[x]),
      unlist(testRanges[x]),
      fragmentLength=NULL,readlength=36
      )
  )
  consensusRanges <- runConsensusRegions(
    testRanges,
    method="majority",
    overlap="any"
    )
  if(unlist(NonPrimaryPeaks["withinsample"])=="drop"){
    consensensusAns <- bplapply(ans,function(x)
      dropNonPrimary(x,consensusRanges)
    )
    ansSummits <- do.call(cbind,bplapply(ans,function(x)
      extractSummits(x,consensusRanges)
    ))
    ansSummitScores <- do.call(cbind,bplapply(ans,function(x)
      extractScores(x,consensusRanges)
    )) 
    if(unlist(NonPrimaryPeaks["betweensample"])=="mean"){
      meanSummits <- rowMeans(ansSummits,na.rm=TRUE)
    }
    if(unlist(NonPrimaryPeaks["betweensample"])=="weightedmean"){
      #meanSummits <- rowMeans(ansSummits,na.omit=TRUE)
      meanSummits <- round(sapply(seq(1,nrow(ansSummits)),function(x)weighted.mean(ansSummits[x,],ansSummitScores[x,])))
    }
  }  
  start(consensusRanges) <- end(consensusRanges) <- meanSummits
  return(consensusRanges)  
}

extractSummits <- function(x,consensusRanges){
  summits <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  summits[overMat[,1]] <- as.vector(start(x[overMat[,2]]))
  return(summits)
}
extractScores <- function(x,consensusRanges,score="summitScores"){
  scores <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  scores[overMat[,1]] <- as.vector(elementMetadata(x[overMat[,2]])[,"summitScores"])
  return(scores)
}

dropNonPrimary <- function(x,consensusRanges,id="elementMetadata.ID",score="summitScores"){
  mat <- as.matrix(findOverlaps(consensusRanges,x))
  tempConsensusRanges <- consensusRanges[mat[,1],]
  elementMetadata(tempConsensusRanges) <- elementMetadata(x[mat[,2]])[,c(id,score)]
  tempConsensusRanges <- tempConsensusRanges[order(elementMetadata(tempConsensusRanges)[,score],decreasing=T),]  
  primaryIDs <- elementMetadata(tempConsensusRanges[match(unique(tempConsensusRanges[,id]),tempConsensusRanges[,id])])[,id]
  x <- x[elementMetadata(x)[,id] %in% primaryIDs]
  return(x)
}

summitPipeline <- function(reads,peakfile,fragmentLength,readlength){
  message("Reading in peaks..",appendLF=FALSE)
  testRanges <- ChIPQC:::GetGRanges(peakfile)
  message("done")  
  if(class(reads) == "GAlignments"){
    message("Alignments loaded")
    ChrLengths <- seqlengths(reads)
  }
  if(class(reads) == "character"){
    message("Reading in alignments..",appendLF=FALSE)
    ChrLengths <- scanBamHeader(reads)[[1]]$targets
    reads <- readGAlignmentsFromBam(reads)
    message("done") 
  }
  message("Calculating fragmentlength..",appendLF=FALSE)
  ccscores <- getShifts(reads,ChrLengths,shiftWindowStart=1,shiftWindowEnd=400)
  fragmentLength <- getFragmentLength(ccscores,readlength)
  message("done")
  message("Extending reads..",appendLF=FALSE)  
  reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
  message("done")
  message("Finding summit locations..",appendLF=FALSE)
  peaks <- runFindSummit(testRanges,reads,fragmentLength=NULL)
  message("done")
  message("Scoring summits..",appendLF=FALSE)
  peaks <- getSummitScore(reads,peaks,fragmentLength=NULL,score="height")
  message("done")
  return(peaks)
}

runConsensusRegions <- function(testRanges,method="majority",overlap="any"){
    if(class(testRanges) == "GRangesList" & length(testRanges) > 1){
      
      reduced <- reduce(unlist(testRanges))
      consensusIDs <- paste0("consensus_",seq(1,length(reduced)))
      elementMetadata(reduced) <- 
      do.call(cbind,lapply(testRanges,function(x)(reduced %over% x)+0))
      if(method=="majority"){
        reducedConsensus <- reduced[rowSums(as.data.frame(elementMetadata(reduced))) > length(testRanges)/2,]
      }
      if(method=="none"){
        reducedConsensus <- reduced
      }
    consensusIDs <- paste0("consensus_",seq(1,length(reducedConsensus)))
    elementMetadata(reducedConsensus) <- cbind(as.data.frame(elementMetadata(reducedConsensus)),consensusIDs)
    return(reducedConsensus)
    
  }
}
runGetShifts <- function(reads,ChrLengths,ChrOfInterestshift,shiftWindowStart=1,shiftWindowEnd=400){

reads <- reads
ChrLengths <- seqlengths(reads)
PosCoverage <- coverage(IRanges(start(reads[strand(reads)=="+"]),start(reads[strand(reads)=="+"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
NegCoverage <- coverage(IRanges(end(reads[strand(reads)=="-"]),end(reads[strand(reads)=="-"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
message("Calculating shift for ",ChrOfInterestshift,"\n")
ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,RleSumAny, verbose = TRUE)         
return(ShiftsTemp)
}
getShifts <- function(reads,ChrLengths,
                      shiftWindowStart=1,shiftWindowEnd=400){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }  
shiftMat <- do.call(cbind,bplapply(names(ChrLengths),function(x)
runGetShifts(reads[seqnames(reads) %in% x],ChrLengths,x,
          shiftWindowStart=1,shiftWindowEnd=400)))
cc_scores <- (rowSums(shiftMat)[1]-rowSums(shiftMat))/rowSums(shiftMat)[1]
return(cc_scores)
}

getFragmentLength <- function(x,readLength){
  #peaks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0)+1
  MaxShift <- which.max(caTools:::runmean(x[-seq(1,(2*readLength))],10))+2*readLength
  
}
  
runFindSummit <- function(testRanges,reads,fragmentLength=NULL){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=F)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
    unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(testRanges))],
    function(x) 
    ChIPQC:::findCovMaxPos(reads[seqnames(testRanges) %in% x],testRanges[seqnames(testRanges) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x],fragmentLength)
    )
  )
  return(test)                                        
}

getSummitScore <- function(reads,summits,fragmentLength=NULL,score="height"){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=F)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
                    as.vector(unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(summits))]),
                    function(x) 
                      runGetSummitScore(reads[seqnames(reads) %in% x],summits[seqnames(summits) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x])
                  )
  )
  return(test)                                         
}
runGetSummitScore <- function(reads,summits,ChrOfInterestshift,FragmentLength=150,score="height"){
    

  cov <- coverage(reads)
    if(score=="height"){
      summitScores <- as.vector(unlist(cov[summits],use.names=F))
    }
  elementMetadata(summits) <- cbind(as.data.frame(elementMetadata(summits)),summitScores)
    return(summits)
}

RleSumAny <- function (e1, e2)
{
  len <- length(e1)
  stopifnot(len == length(e2))
  x1 <- runValue(e1); s1 <- cumsum(runLength(e1))
  x2 <- runValue(e2); s2 <- cumsum(runLength(e2))
  .Call("rle_sum_any",
        as.integer(x1), as.integer(s1),
        as.integer(x2), as.integer(s2),
        as.integer(len),
        PACKAGE = "chipseq")
}
library(BiocParallel)
library(GenomicAlignments)


library(ChIPQC)
36
reads <- "/Users/tcarroll/Downloads/CTCF.cycling.th1DupMarked.bam"
reads <- "/Users/tcarroll/Downloads/DP_CTCFDupMarked.bam"
peakfile <- "/Users/tcarroll/Downloads/DP_CTCF_WithInput_DP_Input_peaks.bed"
testRanges <- ChIPQC:::GetGRanges(peakfile)
ChrLengths <- scanBamHeader(reads)[[1]]$targets
ccscores <- getShifts(reads,ChrLengths,shiftWindowStart=1,shiftWindowEnd=400)
fragmentLength <- getFragmentLength(ccscores,36)
peaks <- runFindSummit(testRanges,reads,fragmentLength)
plott <- regionPlot(reads,peaks,fragmentLength)
plott2 <- regionPlot(reads,testRanges,fragmentLength)
plot(colMeans(plott)[1400:1600],type="l")
lines(colMeans(plott2)[1400:1600],col="red")
bamFiles <- list("/Users/tcarroll/Downloads/CTCFCyclingTHDupMarked.bam","/Users/tcarroll/Downloads/DP_CTCFDupMarked.bam")
testRanges <- list("/Users/tcarroll/Downloads/CTCFCyclingTH_WithInput_InputCyclingTH_peaks.bed","/Users/tcarroll/Downloads/DP_CTCF_WithInput_DP_Input_peaks.bed")
names(testRanges) <- c("Cycling","DP")
names(bamFiles) <- c("Cycling","DP") 
consensusRegions <- findconsensusRegions(testRanges,bamFiles)
tt2 <- ChIPQC:::GetGRanges(testRanges[[2]])[ChIPQC:::GetGRanges(testRanges[[2]]) %over% consensusRegions]
tt1 <- ChIPQC:::GetGRanges(testRanges[[1]])[ChIPQC:::GetGRanges(testRanges[[1]]) %over% consensusRegions]
test <- bplapply(bamFiles,
            function(x)
            regionPlot(x,consensusRegions,style="point",format="bam",FragmentLength=130)
)

test2 <- bplapply(bamFiles,
                  function(x)
                    regionPlot(x,ChIPQC:::GetGRanges(testRanges[[2]]),style="point",format="bam",FragmentLength=130)
)
test22 <- bplapply(bamFiles,
                  function(x)
                    regionPlot(x,ChIPQC:::GetGRanges(testRanges[[1]]),style="point",format="bam",FragmentLength=130)
)

testset <- reduce(c(ChIPQC:::GetGRanges(testRanges[[2]]),ChIPQC:::GetGRanges(testRanges[[1]])))

test3 <- bplapply(bamFiles,
                  function(x)
                    regionPlot(x,testset,style="point",format="bam",FragmentLength=130)
)

testtt1 <- bplapply(bamFiles,
                  function(x)
                    regionPlot(x,tt1,style="point",format="bam",FragmentLength=130)
)
testtt2 <- bplapply(bamFiles,
                    function(x)
                      regionPlot(x,tt2,style="point",format="bam",FragmentLength=130)
)

### Is the ordering creating a skew again for average plot in selecting summit!
par(mfrow=c(3,2))
plot(colMeans(unlist(test[[1]])),type="l")
plot(colMeans(unlist(test[[2]])),type="l")
plot(colMeans(unlist(test2[[1]])),type="l")
plot(colMeans(unlist(test2[[2]])),type="l")
plot(colMeans(unlist(test3[[1]])),type="l")
plot(colMeans(unlist(test3[[2]])),type="l")
#testset <- reduce(c(ChIPQC:::GetGRanges(testRanges[[2]]),ChIPQC:::GetGRanges(testRanges[[1]])))
dev.off()
plot(colMeans(unlist(test[[1]])),type="l")
lines(colMeans(unlist(test3[[1]])),col="red")
lines(colMeans(unlist(test2[[1]])),col="blue")
lines(colMeans(unlist(test22[[1]])),col="blue",lty=2)

plot(colMeans(unlist(test[[2]])),type="l")
lines(colMeans(unlist(test3[[2]])),col="red")
lines(colMeans(unlist(test2[[2]])),col="blue")
lines(colMeans(unlist(test22[[2]])),col="blue",lty=2)
lines(colMeans(unlist(testtt2[[2]])),col="darkgreen",lty=2)
lines(colMeans(unlist(testtt1[[2]])),col="purple",lty=2)

plot(colMeans(unlist(test[[1]]))[1480:1520],type="l")
lines(colMeans(unlist(test3[[1]]))[1480:1520],col="red")
lines(colMeans(unlist(test2[[1]])),col="blue")
lines(colMeans(unlist(test22[[1]])),col="blue",lty=2)
lines(colMeans(unlist(testtt2[[1]])),col="darkgreen",lty=2)
lines(colMeans(unlist(testtt1[[1]])),col="purple",lty=2)
