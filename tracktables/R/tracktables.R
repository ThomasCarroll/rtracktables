#' Plot coverage of points or regions.
#'
#' @param bamFile Character vector for location of BAM file.
#' @param testRanges GRanges object of regions to plot.
#' @param FragmentLength Integer vector Predicted or expected fragment length.
#' @param style Point or region (see details)
#' @param distanceAround Distance around centre of region to be used for plotting
#' @param distanceInRegionStart Distance into region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionStart Distance out from region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceInRegionEnd Distance into region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionEnd Distance out from region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param paired Is data paired end 
#' @param normalize Calculate coverage as RPM. Presently only RPM available.
#' @param plotBy Score to be used for plotting. Presently only coverage.
#' @param removeDup Remove duplicates before calculating coverage.
#' @param format BAM or BigWig
#' @param seqlengths Chromosomes to be used. If missing will report all.
#' @param forceFragment Centre fragment and force consistent fragment width.
#' @return ChIPprofile A ChIPprofile object. 
#' @export
#' @import IRanges GenomicRanges ggplot2 QuasR rtracklayer GenomicAlignments GenomicRanges XVector Rsamtools
#' 
regionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,verbose=T,format="bam",seqlengths=NULL,forceFragment=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  result <- runRegionPlot(bamFile,testRanges,nOfWindows,FragmentLength,style,distanceAround,distanceInRegionStart,distanceOutRegionStart,distanceInRegionEnd,distanceOutRegionEnd,paired,normalize,plotBy,removeDup,format,seqlengths,forceFragment)
  return(result)  
}

runRegionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,format="bam",seqlengths=NULL,forceFragment=NULL){

  #bamFile <- "/home//pgellert/Dropbox (Lymphocyte_Developme)/WeiWeiLiang/RNAPII/Sample_R1-0hDupMarked.bam"
  #bamFile <-"Downloads//mergedETOH.bwRange5.bw"
  #bamFile <-"Downloads//Sample_R1-6hDupMarkedNormalised.bw"
  #testRanges <- mm9PC
  #nOfWindows=100
  #FragmentLength=150
  #style="region"
  #distanceAround=1500
  #distanceInRegionStart=1500
  #distanceOutRegionStart=1500
  #distanceInRegionEnd=1500
  #distanceOutRegionEnd=1500
  #paired=F
  #normalize="RPM"
  #plotBy="coverage"
  #removeDup=F  
  #format="bigwig"
  #seqlengths=NULL
  
  ## Initialize empty matricies and paramaters for collecting coverage analysis
  if(style == "region" | style=="regionandpoint"){
    posRegionStartMat <- NULL
    posRegionEndMat <- NULL
    negRegionStartMat <- NULL
    negRegionEndMat <- NULL
    RegionsMat <- NULL
    maxDistance <- max(distanceOutRegionStart,distanceOutRegionEnd)
    distanceUpStart <- distanceOutRegionStart
    distanceDownEnd <- distanceOutRegionEnd
    
  }
  
  if(style == "point"){
    PosRegionMat <- NULL
    NegRegionMat <- NULL
    RegionsMat <- NULL    
    maxDistance=distanceAround
    distanceUpStart <- distanceAround
    distanceDownEnd <- distanceAround    
  }
  
  if(format == "bam"){
    ## Get all chromosomes in bamFile
    message("Reading Bam header information...",appendLF = FALSE)
    allchrs <- names(scanBamHeader(bamFile)[[1]]$targets)
    lengths <- as.vector(scanBamHeader(bamFile)[[1]]$targets)
    names(lengths) <- allchrs
    message("..Done")
  }
  if(format=="bigwig"){
    message("Importing BigWig...",appendLF = FALSE)
    genomeCov <- import.bw(bamFile,as = "RleList")
    if(is.null(seqlengths)){
      seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
    }else{
      seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }
  
  
  ## Filter testRanges to those contained within chromosomes.
  message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
  testRangeNames <- unique(seqnames(testRanges))
  temptestranges <- GRanges()
  for(i in 1:length(testRangeNames)){
    perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
    temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                    & start(perchrRanges)-maxDistance > 0 ])
    #print(i)
  }
  
  message("..Done")
  message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
  testRanges <- temptestranges
  message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
  
  strand(testRanges[strand(testRanges) == "*"]) <- "+"
  testRangesPos <- testRanges[strand(testRanges) == "+"]
  testRangesNeg <- testRanges[strand(testRanges) == "-"]
  message("..Done")
  if(style=="region"){
    message("Filtering regions which are smaller than windows into region...",appendLF = FALSE)
    ## Split Regions into those on positive and negative strands..
    testRangesPos <- testRangesPos[(end(testRangesPos)-distanceInRegionEnd) - (start(testRangesPos)+distanceInRegionStart) > nOfWindows]
    testRangesNeg <- testRangesNeg[(end(testRangesNeg)-distanceInRegionStart) - (start(testRangesNeg)+distanceInRegionEnd) > nOfWindows]
    message("..Done")
  }  
  
  message("Found ",length(testRangesPos)," Watson strand regions")
  message("Found ",length(testRangesNeg)," Crick strand regions")
  
  
  ## Extend regions and get positive versus negative reads.
  message("Extending regions..",appendLF=FALSE)    
  exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStart,end(testRangesPos)+distanceDownEnd)),
                     GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEnd,end(testRangesNeg)+distanceUpStart))
  )
  message("...done")   
  
  
  if(!removeDup){
    Param <- ScanBamParam(which=GRanges(seqnames=seqnames(exttestRanges[seqnames(exttestRanges) %in% allchrs]),IRanges(start=start(exttestRanges[seqnames(exttestRanges) %in% allchrs]),end=end(exttestRanges[seqnames(exttestRanges) %in% allchrs]))))
  }else{
    Param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),which=GRanges(seqnames=seqnames(exttestRanges[seqnames(exttestRanges) %in% allchrs]),IRanges(start=start(exttestRanges[seqnames(exttestRanges) %in% allchrs]),end=end(exttestRanges[seqnames(exttestRanges) %in% allchrs]))))
  }
  
  if(format == "bam"){
    message("Reading tags from ",bamFile,appendLF=FALSE)
    totalReads <- alignmentStats(bamFile)[,"mapped"]
    if(paired==FALSE){
      total <- readGAlignmentsFromBam(bamFile,param=Param)
      message("..Done.\nRead in ",length(total)," reads")
      message("Extending reads to fragmentlength of ",FragmentLength,appendLF=F)
      temp <- resize(as(total,"GRanges"),FragmentLength,"start")
      message("..done")
    }
    if(paired==TRUE){
      #bamFile <- BamFile(bamFile,yieldSize=50000000L)
      #tempPaired <- readGAlignmentPairsFromBam(bamFile,param=Param)
      #tempPaired <- tempPaired[isProperPair(tempPaired)]
      gaPaired <- readGAlignmentsFromBam(bamFile, 
                                         param=ScanBamParam(what=c("mpos"),
                                                            flag=scanBamFlag(isProperPair = TRUE,isFirstMateRead = TRUE)))      
      tempPos <- GRanges(seqnames(gaPaired[strand(gaPaired) == "+"]),
                         IRanges(
                           start=start(gaPaired[strand(gaPaired) == "+"]),
                           end=elementMetadata(gaPaired[strand(gaPaired) == "+"])$mpos
                           +qwidth(gaPaired[strand(gaPaired) == "+"])))
      tempNeg <- GRanges(seqnames(gaPaired[strand(gaPaired) == "-"]),
                         IRanges(
                           start=elementMetadata(gaPaired[strand(gaPaired) == "-"])$mpos,                        
                           end=end(gaPaired[strand(gaPaired) == "-"])
                         )) 
      temp <- c(tempPos,tempNeg)                
      #temp <- GRanges(seqnames(tempPaired),IRanges(start(left(tempPaired)),end(right(tempPaired))))
      message("..Done.\nRead in ",length(temp)," reads")
      if(!is.null(forceFragment)){
        message("Forcing fragments to be centred and set to ",forceFragment,"..",appendLF=FALSE)
        temp <- resize(temp,forceFragment,"center")
        message("..done")        
      }
      message("..done")
    }  
    message("Calculating coverage..",appendLF=FALSE)
    genomeCov <- coverage(temp)
    message("..done")
  }
  chromosomes <- seqlevels(genomeCov) 
  
  if(style=="point"){
    testRangesPos <- resize(testRangesPos,1,"center")
    testRangesNeg <- resize(testRangesNeg,1,"center")
    RangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStart,start(testRangesPos)+distanceDownEnd),name=testRangesPos$name)
    RangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(end(testRangesNeg)-distanceDownEnd,end(testRangesNeg)+distanceUpStart),name=testRangesNeg$name)  
    
    for(c in 1:length(chromosomes)){
      if(length(RangesPos[seqnames(RangesPos) %in% chromosomes[c]]) > 0){
        PosRegionMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesPos[seqnames(RangesPos) %in% chromosomes[c]])]),ncol=mean(width(RangesPos)),byrow=TRUE)
        rownames(PosRegionMat) <- RangesPos[seqnames(RangesPos) %in% chromosomes[c]]$name
        #      print("done1.3")
      }
      if(length(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]) > 0){
        NegRegionMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]])])),ncol=mean(width(RangesNeg)),byrow=TRUE)
        rownames(NegRegionMat) <- RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]$name
      }    
      RegionsMat <- rbind(RegionsMat,PosRegionMat,NegRegionMat)
    }
    return(RegionsMat)
  }
  if(style=="region"){
    
    message("Defining flanks of regions..",appendLF=FALSE)
    
    startRegionRangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceOutRegionStart,start(testRangesPos)+distanceInRegionStart),name=testRangesPos$name)
    endRegionRangesPos <- GRanges(seqnames(testRangesPos),IRanges(end(testRangesPos)-distanceInRegionEnd,end(testRangesPos)+distanceOutRegionEnd),name=testRangesPos$name)
    startRegionRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(end(testRangesNeg)-distanceInRegionStart,end(testRangesNeg)+distanceOutRegionStart),name=testRangesNeg$name)
    endRegionRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceOutRegionEnd,start(testRangesNeg)+distanceInRegionEnd),name=testRangesNeg$name)
    
    testRangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)+distanceInRegionStart,end(testRangesPos)-distanceInRegionEnd),name=testRangesPos$name)
    testRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)+distanceInRegionEnd,end(testRangesNeg)-distanceInRegionStart),name=testRangesNeg$name)     
    message("...Done")
    meansList <- vector("numeric")
    grListWindowsPos <- GRanges()
    grListWindowsNeg <- GRanges()
    #grListWindows <- list()
    message("Making windows..",appendLF=FALSE)
    
    
    
    if(length(testRangesNeg) > 0){
      grWidths <- width(testRangesPos)  
      windows <- floor(grWidths%/%nOfWindows)
      extraLastWindow <- grWidths%%nOfWindows
      addToWindow <- 0
      for(i in 1:nOfWindows){
        
        if(i == nOfWindows){
          addToWindow <- extraLastWindow 
        }
        grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
          start(testRangesPos)+(windows*(i-1)),
          start(testRangesPos)+(windows*i)-1+addToWindow),name=testRangesPos$name))
      }
    }
    grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$name)]
    if(length(testRangesNeg) > 0){
      grWidths <- width(testRangesNeg)
      windows <- floor(grWidths%/%nOfWindows)
      extraLastWindow <- grWidths%%nOfWindows
      addToWindow <- 0
      for(i in 1:nOfWindows){
        
        if(i == nOfWindows){
          addToWindow <- extraLastWindow 
        }
        grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(     
          end(testRangesNeg)-(windows*i)+1-addToWindow,
          end(testRangesNeg)-(windows*(i-1))),name=testRangesNeg$name))
      }
      grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$name)]
    }
    grListWindows <- c(grListWindowsPos,grListWindowsNeg)
    message("..done\n")
    message("Calculating per chromsome")      
    for(c in 1:length(chromosomes)){
      
      message("Processing inner region windows in ",chromosomes[c])
      covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
      doubleTemp <- viewMeans(covPerPeak)
      names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$name)
      meansList <- c(meansList,doubleTemp)
      message("..done")
      message("Processing flanking windows in ",chromosomes[c])      
      
      tempstartRegionRangesPosMat <- NULL
      tempendRegionRangesPosMat <- NULL  
      tempstartRegionRangesNegMat <- NULL
      tempendRegionRangesNegMat <- NULL    
      
      if(length(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]) > 0){
        tempstartRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]])]),ncol=mean(width(startRegionRangesPos)),byrow=TRUE)
        rownames(tempstartRegionRangesPosMat) <- startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]$name
      }
      
      if(length(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]) > 0){
        tempendRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]])]),ncol=mean(width(endRegionRangesPos)),byrow=TRUE)
        rownames(tempendRegionRangesPosMat) <- endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]$name
      }
      if(length(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]) > 0){
        tempstartRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]])])),ncol=mean(width(startRegionRangesNeg)),byrow=TRUE)
        rownames(tempstartRegionRangesNegMat) <- rev(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]$name)
      }
      #print("done4")  
      if(length(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]) > 0){
        #tempNegTTSMat <- matrix(rev(as.vector(unlist(Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(ttsRangesNeg[seqnames(ttsRangesNeg) %in% chromosomes[c]]))))),ncol=mean(width(ttsRangesNeg)),byrow=TRUE)
        tempendRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]])])),ncol=mean(width(endRegionRangesNeg)),byrow=TRUE)
        rownames(tempendRegionRangesNegMat) <- rev(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]$name)
      }
      #print("done5")
      
      posRegionStartMat <- rbind(posRegionStartMat,tempstartRegionRangesPosMat)
      posRegionEndMat <- rbind(posRegionEndMat,tempendRegionRangesPosMat)
      negRegionStartMat <- rbind(negRegionStartMat,tempstartRegionRangesNegMat)
      negRegionEndMat <- rbind(negRegionEndMat,tempendRegionRangesNegMat)
      tempstartRegionRangesPosMat <- NULL
      tempendRegionRangesPosMat <- NULL  
      tempstartRegionRangesNegMat <- NULL
      tempendRegionRangesNegMat <- NULL           
      message("..done")
    }
    AllRegionStart <- rbind(posRegionStartMat,negRegionStartMat)
    AllRegionEnd <- rbind(posRegionEndMat,negRegionEndMat)
    meansMat <- matrix(meansList,ncol=nOfWindows,byrow=T)
    rownames(meansMat) <- matrix(names(meansList),ncol=nOfWindows,byrow=T)[,1]
    start <- cbind(seq(1:length(colMeans(AllRegionStart))),colMeans(AllRegionStart))
    mid <- cbind(max(start[,1])+seq(1:length(colMeans(meansMat)))*100,colMeans(meansMat))
    end <- cbind(max(mid[,1])+seq(1:length(colMeans(AllRegionEnd))),colMeans(AllRegionEnd))
    #reportRanges <- testRanges[match(rownames(meansMat),testRanges$name)]
    #returmm9PC[match(mm9PC$name,rownames(temp))]n(list(meansMat,AllRegionStart,AllRegionEnd,rbind(start,mid,end)))
    #reportRanges <- testRanges[match(rownames(meansMat),testRanges$name)]
    #elementMetadata(reportRanges) <- cbind(elementMetadata(reportRanges),AllRegionStart,meansMat,AllRegionEnd)
    #new("ChIPprofile",reportRanges,profile=list())
    #return(profiles)
    profileMat <- cbind(AllRegionStart,meansMat,AllRegionEnd)
    colnames(profileMat) <- c(paste0("Region_Start",seq(0-distanceOutRegionStart,-1)),"Region_Start",paste0("Region_Start",seq(1,distanceOutRegionStart)),
    paste0(seq(1:100),"%_ofRegion"),
    paste0("Region_End",seq(0-distanceOutRegionEnd,-1)),"Region_End",paste0("Region_End",seq(1,distanceOutRegionEnd))
    )
    filteredRanges <- c(testRangesPos,testRangesNeg)
    profileSample <- SummarizedExperiment(profileMat,rowData=filteredRanges[match(rownames(profileMat),filteredRanges$name)])
    return(new("ChIPprofile",profileSample,profile=list()))
  } 
}

