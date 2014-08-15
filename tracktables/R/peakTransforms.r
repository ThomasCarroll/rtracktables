findconsensusRegions <- function(bamFile,testRanges,method="majority",summit="weightedmean",resizepeak="asw",overlap="any",fragmentLength=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  consensusRanges <- runConsensusRegions(bamFile,testRanges,method="majority",summit="weightedmean",resizepeak="asw",overlap="any",fragmentLength=NULL)
  
  return(result)  
}

runConsensusRegions(testRanges,method="majority"){
  if(class(testRanges) == "GRangesList" & length(testRanges) > 1){
    
    reduced <- reduce(unlist(GRangesList(testRanges)))
    elementMetadata(reduced) <- 
    do.call(cbind,lapply(testRanges,function(x)(reduced %over% x)+0))
    if(method="majority"){
      reducedConsensus <- reduced[rowSums(as.data.frame(elementMetadata(reduced))) > length(testRanges)/2,]
    }
    if(method="none"){
      reducedConsensus <- reduced
    }
    
  }
  return(reducedConsensus)
  
}

runGetShifts <- function(reads,ChrLengths,ChrOfInterestshift,WindowStart,shiftWindowStart=1,shiftWindowEnd=400){
reads <- reads[seqnames(reads) %in% ChrOfInterestshift]
ChrLengths <- seqlengths(reads)
PosCoverage <- coverage(IRanges(start(reads[strand(reads)=="+"]),start(reads[strand(reads)=="+"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
NegCoverage <- coverage(IRanges(end(reads[strand(reads)=="-"]),end(reads[strand(reads)=="-"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
message("Calculating shift for ",ChrOfInterestshift,"\n")
ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,RleSumAny, verbose = TRUE)         
return(ShiftsTemp)
}
getShifts <- function(reads,ChrLengths,x,WindowStart,
                      shiftWindowStart=1,shiftWindowEnd=400){
shiftMat <- do.call(cbind,bplapply(names(ChrLengths),function(x)
runGetShifts(reads,ChrLengths,x,WindowStart,
          shiftWindowStart=1,shiftWindowEnd=400)))
cc_scores <- (rowSums(shiftMat)[1]-rowSums(shiftMat))/rowSums(shiftMat)[1])
}

getShifts <- function(reads,ChrLengths,x,WindowStart,
                      shiftWindowStart=1,shiftWindowEnd=400)
  
funFindSummit <- function(testRanges,bamFile,FragmentLength){
  total <- readGAlignmentsFromBam(bamFile)
  message("..Done.\nRead in ",length(total)," reads")
  message("Extending reads to fragmentlength of ",FragmentLength,appendLF=F)
  temp <- resize(as(total,"GRanges"),FragmentLength,"start")
  rm(total)
  test <- do.call(c,
                  bplapply(
    unique(seqnames(temp))[unique(seqnames(temp)) %in% unique(seqnames(testRanges))],
    function(x) 
    ChIPQC:::findCovMaxPos(temp[seqnames(testRanges) %in% x],testRanges[seqnames(testRanges) %in% x],seqlengths(temp)[names(seqlengths(temp)) %in% x],FragmentLength)
    )
  )
                                          
}