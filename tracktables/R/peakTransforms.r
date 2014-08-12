consensusRegions <- function(bamFile,testRanges,method="majority",summit="weightedmean",resizepeak="asw",overlap="any",fragmentLength=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  consensusRanges <- findConsensusRegions(bamFile,testRanges,method="majority",summit="weightedmean",resizepeak="asw",overlap="any",fragmentLength=NULL)
  
  return(result)  
}

findConsensusRegions(testRanges,method="majority"){
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
findSummit <- function(testRanges,bamFile,FragmentLength){
  total <- readGAlignmentsFromBam(bamFile,param=Param)
  message("..Done.\nRead in ",length(total)," reads")
  message("Extending reads to fragmentlength of ",FragmentLength,appendLF=F)
  temp <- resize(as(total,"GRanges"),FragmentLength,"start")
  rm(total)
  cov <- coverage(temp)
  rm(total)
}