
GetAroundTSS <- function(FullGeneBounds,distance,distanceIn){
  require(GenomicRanges)
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  TempFramePos <- cbind(as.data.frame(elementMetadata(PosGenes)),width(PosGenes))
  TempFrameNeg <- cbind(as.data.frame(elementMetadata(NegGenes)),width(NegGenes))
  colnames(TempFramePos)[3] <- "WidthOfOriginalGene"
  colnames(TempFrameNeg)[3] <- "WidthOfOriginalGene"  
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distance,end=(as.vector(start(ranges(PosGenes))))+distanceIn),strand=strand(PosGenes),TempFramePos)
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes)))-distanceIn),end=(as.vector(end(ranges(NegGenes))))+distance),strand=strand(NegGenes),TempFrameNeg)
  names(NewPosGenes) <- names(PosGenes)
  names(NewNegGenes) <- names(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  return(AllPromoters)
}

GetTargetTSS <- function(RegionRanges,distance,distanceIn=0,min){
  RegionRanges <- RegionRanges[width(RegionRanges) > min]
  TSSs <- GetAroundTSS(RegionRanges,distance,distanceIn=distanceIn)
  #print(TSSs)
  #  Centre <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)-distance),end=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)+distance)),strand="*")
  #RegionRangesList <- GRangesList(TSSs[strand(TSSs) == "+"],TSSs[strand(TSSs) == "-"])
  return(TSSs)
}

GetAroundTTS <- function(FullGeneBounds,distance,distanceIn){
  require(GenomicRanges)
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  TempFramePos <- cbind(as.data.frame(elementMetadata(PosGenes)),width(PosGenes))
  TempFrameNeg <- cbind(as.data.frame(elementMetadata(NegGenes)),width(NegGenes))
  colnames(TempFramePos)[3] <- "WidthOfOriginalGene"
  colnames(TempFrameNeg)[3] <- "WidthOfOriginalGene"  
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(end(ranges(PosGenes))))-distanceIn,end=(as.vector(end(ranges(PosGenes))))+distance),strand=strand(PosGenes),TempFramePos)
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(start(ranges(NegGenes)))-distance),end=(as.vector(start(ranges(NegGenes))))+distanceIn),strand=strand(NegGenes),TempFrameNeg)
  names(NewPosGenes) <- names(PosGenes)
  names(NewNegGenes) <- names(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  return(AllPromoters)
}

GetTargetTTS <- function(RegionRanges,distance,distanceIn=0,min){
  RegionRanges <- RegionRanges[width(RegionRanges) > min]
  TTSs <- GetAroundTTS(RegionRanges,distance,distanceIn=distanceIn)
  #print(TSSs)
  #  Centre <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)-distance),end=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)+distance)),strand="*")
  #RegionRangesList <- GRangesList(TSSs[strand(TSSs) == "+"],TSSs[strand(TSSs) == "-"])
  return(TTSs)
}


GetGene <- function(FullGeneBounds,distanceFromStart,distanceFromEnd,min){
  require(GenomicRanges)
  FullGeneBounds <- FullGeneBounds[width(FullGeneBounds) > min]
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distanceFromStart,end=(as.vector(end(ranges(PosGenes))))+distanceFromEnd),strand=strand(PosGenes),elementMetadata(PosGenes))
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(start(ranges(NegGenes)))-distanceFromEnd),end=(as.vector(end(ranges(NegGenes))))+distanceFromStart),strand=strand(NegGenes),elementMetadata(NegGenes))
  names(NewPosGenes) <- names(PosGenes)
  names(NewNegGenes) <- names(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  return(AllPromoters)
}
regionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,verbose=T,format="bam",seqlengths=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  result <- runRegionPlot(bamFile,testRanges,nOfWindows,FragmentLength,style,distanceAround,distanceInRegionStart,distanceOutRegionStart,distanceInRegionEnd,distanceOutRegionEnd,paired,normalize,plotBy,removeDup,format,seqlengths)
  return(result)  
}
runRegionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,format="bam",seqlengths=NULL){
  require(QuasR)
  require(rtracklayer)  
  require(GenomicAlignments)
  require(GenomicRanges)  
  
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
    perchrRanges <- testRanges[seqnames(testRanges) %in% testRangeNames[i]]
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
    tempPaired <- readGAlignmentPairsFromBam(bamFile)
    tempPaired <- tempPaired[isProperPair(tempPaired)]
    temp <- GRanges(seqnames(tempPaired),IRanges(start(left(tempPaired)),end(right(tempPaired))))
    message("..Done.\nRead in ",length(temp)," reads")
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
    return(list(meansMat,AllRegionStart,AllRegionEnd,rbind(start,mid,end)))

} 
}

mm9Genes <- read.delim("/Users/tcarroll/Downloads/mm9Genes_May2012.txt",sep="\t",h=T)

mm9GeneRanges <- GRanges(seqnames=paste0("chr",mm9Genes[,3]),ranges=IRanges(start=mm9Genes[,1],end=mm9Genes[,2]),strand=mm9Genes[,4],name=mm9Genes[,5],biotype=mm9Genes[,6])
JustChrOfInterest <- unique(as.vector(seqnames(mm9GeneRanges)))[grep("chr\\d.|chr\\d|chrX|chrY|chrM",unique(as.vector(seqnames(mm9GeneRanges))))]
mm9PC <- mm9GeneRanges[mm9GeneRanges$biotype == "protein_coding"]
mm9PC <- mm9PC[order(width(mm9PC),decreasing=T)]
mm9PC <- mm9PC[match(unique(mm9PC$name),mm9PC$name)]
mm9PC <- mm9PC[!mm9PC$name == ""]
mm9PC <- mm9PC[seqnames(mm9PC) %in% JustChrOfInterest]
mm9PC <- mm9PC[width(mm9PC) > 200]

#mm9DF <- as.data.frame(mm9PC)

mm9TSS <- GetTargetTSS(mm9PC,300,300,500)
#bamFile <-"Downloads//mergedETOH.bwRange5.bw"
#bamFile <-"Downloads//Sample_R1-6hDupMarkedNormalised.bw"
#res <- regionPlot("/home//pgellert/Dropbox (Lymphocyte_Developme)/WeiWeiLiang/RNAPII/Sample_R1-0hDupMarked.bam",mm9TSS)
#res1 <- regionPlot("/home//pgellert/Dropbox (Lymphocyte_Developme)/WeiWeiLiang/nucelosome/fromMRC/mergedETOH.bwRange5.bw",mm9PC,style="region",format="bigwig")
nuc0h <- regionPlot("Downloads//mergedETOH.bwRange5.bw",mm9PC,style="region",format="bigwig")
nuc6h <- regionPlot("Downloads//mergedOHT.bwRange5.bw",mm9PC,style="region",format="bigwig")
pol6h <- regionPlot("Downloads//Sample_R2-6hDupMarkedNormalised.bw",mm9PC,style="region",format="bigwig")
pol0h <- regionPlot("Downloads//Sample_R1-0hDupMarkedNormalised.bw",mm9PC,style="region",format="bigwig")
nuc0h4 <- regionPlot("Downloads//mergedETOH.bwRange4.bw",mm9PC,style="region",format="bigwig")

nuc6h4 <- regionPlot("Downloads//mergedOHT.bwRange4.bw",mm9PC,style="region",format="bigwig")
nuc6h3 <- regionPlot("Downloads//mergedOHT.bwRange3.bw",mm9PC,style="region",format="bigwig")
nuc6h2 <- regionPlot("Downloads//mergedOHT.bwRange2.bw",mm9PC,style="region",format="bigwig")
dir("Downloads")
b3ranges <- ChIPQC:::GetGRanges("Downloads//B3Peaks.bed.txt")
Ik6h <- regionPlot("Downloads//mergedOHT.bwRange5.bw",b3ranges,style="point",format="bigwig")
Ik0h <- regionPlot("Downloads//mergedETOH.bwRange5.bw",b3ranges,style="point",format="bigwig")
Ik6h4 <- regionPlot("Downloads//mergedOHT.bwRange4.bw",b3ranges,style="point",format="bigwig")
Ik0h4 <- regionPlot("Downloads//mergedETOH.bwRange4.bw",b3ranges,style="point",format="bigwig")
Ik6h3 <- regionPlot("Downloads//mergedOHT.bwRange3.bw",b3ranges,style="point",format="bigwig")
Ik0h3 <- regionPlot("Downloads//mergedETOH.bwRange3.bw",b3ranges,style="point",format="bigwig")
Ik6h2 <- regionPlot("Downloads//mergedOHT.bwRange2.bw",b3ranges,style="point",format="bigwig")
Ik0h2 <- regionPlot("Downloads//mergedETOH.bwRange2.bw",b3ranges,style="point",format="bigwig")


Ik0h2me2 <- regionPlot("Downloads//Mi2b_IKneg.bw",b3ranges,style="point",format="bigwig")
Ik6h2me2 <- regionPlot("Downloads//Mi2b_IKpos.bw",b3ranges,style="point",format="bigwig")

tss0h2me2 <- regionPlot("Downloads//Mi2b_IKneg.bw",mm9PC,style="region",format="bigwig")
tss6h2me2 <- regionPlot("Downloads//Mi2b_IKpos.bw",mm9PC,style="region",format="bigwig")

par(mfrow=c(2,1))
plot(tss6h2me2[[4]],col="darkblue",type="l")
plot(tss0h2me2[[4]],col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h2me2),col="darkblue",type="l")
plot(log2(colMeans(Ik6h2me2)/colMeans(Ik0h2me2)),col="darkred",type="l")
plot(colMeans(Ik0h2me2),col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h),col="darkblue",type="l")
plot(log2(colMeans(Ik6h)/colMeans(Ik0h)),col="darkred",type="l")
plot(colMeans(Ik0h),col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h),col="darkblue",type="l")
plot(log2(colMeans(Ik6h)/colMeans(Ik0h)),col="darkred",type="l")
plot(colMeans(Ik0h),col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h4),col="darkblue",type="l")
plot(log2(colMeans(Ik6h4)/colMeans(Ik0h4)),col="darkred",type="l")
plot(colMeans(Ik0h4),col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h3),col="darkblue",type="l")
plot(log2(colMeans(Ik6h3)/colMeans(Ik0h3)),col="darkred",type="l")
plot(colMeans(Ik0h3),col="darkred",type="l")
###############
par(mfrow=c(3,1))
plot(colMeans(Ik6h2),col="darkblue",type="l")
plot(log2(colMeans(Ik6h2)/colMeans(Ik0h3)),col="darkred",type="l")
plot(colMeans(Ik0h2),col="darkred",type="l")






plot(pol0h[[4]][1:3000,1],runmean(pol0h[[4]][1:3000,2],1),col="darkred",type="l",ylim=c(0,3.5))
lines(nuc6h[[4]][1:3000,1],runmean(nuc6h[[4]][1:3000,2],1)-2.8,col="darkblue")
par(mfrow=c(2,3))
plot(pol0h[[4]][1:3000,1],runmean(pol0h[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc6h[[4]][1:3000,1],runmean(nuc6h[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc6h4[[4]][1:3000,1],runmean(nuc6h4[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc6h3[[4]][1:3000,1],runmean(nuc6h3[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc6h2[[4]][1:3000,1],runmean(nuc6h2[[4]][1:3000,2],1),col="darkred",type="l")
png("Downloads/Profiles.png",width=1000,height=2000)
par(mfrow=c(6,1))
plot(pol0h[[4]][,1],runmean(pol0h[[4]][,2],1),col="darkred",type="l",lwd=2)
plot(nuc6h[[4]][,1],runmean(nuc6h[[4]][,2],1),col="darkred",type="l",lwd=2)
plot(nuc6h4[[4]][,1],runmean(nuc6h4[[4]][,2],1),col="darkred",type="l",lwd=2)
plot(nuc6h3[[4]][,1],runmean(nuc6h3[[4]][,2],1),col="darkred",type="l",lwd=2)
plot(nuc6h2[[4]][,1],runmean(nuc6h2[[4]][,2],1),col="darkred",type="l",lwd=2)
dev.off()

lines(nuc6h[[4]][1:3000,1],runmean(nuc6h[[4]][1:3000,2],1)-2.8,col="darkblue")


res3 <- readGAlignmentsFromBam("/Users/tcarroll/Downloads//mergedETOH.bamRange5.bam", param=ScanBamParam(what=c("mpos"),flag=scanBamFlag(isProperPair = TRUE,isFirstMateRead = TRUE)))

nuc0hrange5 <- regionPlot("/Users/tcarroll/Downloads//mergedETOH.bamRange5.bam",mm9PC,style="region",format="bam",paired=T,forceFragment=80)
nuc0hrange4 <- regionPlot("/Users/tcarroll/Downloads//mergedETOH.bamRange4.bam",mm9PC,style="region",format="bam",paired=T,forceFragment=80)
nuc0hrange3 <- regionPlot("/Users/tcarroll/Downloads//mergedETOH.bamRange3.bam",mm9PC,style="region",format="bam",paired=T,forceFragment=80)
nuc0hrange2 <- regionPlot("/Users/tcarroll/Downloads//mergedETOH.bamRange2.bam",mm9PC,style="region",format="bam",paired=T,forceFragment=50)

dir2 <- dir("/Users/tcarroll/Downloads//",pattern="mergedETOH.bamRange\\d.bam$",full.names=T)
mclapply(dir2,indexBam)
regionPlot("/Users/tcarroll/Downloads//mergedETOH.bamRange5.bam",
           mm9PC,
           style="region",
           format="bam",paired=T,
           forceFragment=80)
par(mfrow=c(2,2))
plot(nuc0hrange5[[4]][1:3000,1],runmean(nuc0hrange5[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc0hrange4[[4]][1:3000,1],runmean(nuc0hrange4[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc0hrange3[[4]][1:3000,1],runmean(nuc0hrange3[[4]][1:3000,2],1),col="darkred",type="l")
plot(nuc0hrange2[[4]][1:3000,1],runmean(nuc0hrange2[[4]][1:3000,2],1),col="darkred",type="l")

par(mfrow=c(2,2))
plot(nuc0hrange5[[4]][,1],runmean(nuc0hrange5[[4]][,2],1),col="darkred",type="l")
plot(nuc0hrange4[[4]][,1],runmean(nuc0hrange4[[4]][,2],1),col="darkred",type="l")
plot(nuc0hrange3[[4]][,1],runmean(nuc0hrange3[[4]][,2],1),col="darkred",type="l")
plot(nuc0hrange2[[4]][,1],runmean(nuc0hrange2[[4]][,2],1),col="darkred",type="l")

