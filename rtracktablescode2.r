MakeIGVSampleMetadata <- function(sampleMetadata,SampleSheet,igvdirectory){
  write.table("#sampleTable",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  colnames(sampleMetadata)[1] <- "Linking_id"
  print(colnames(sampleMetadata))
  sampleMetadata <- as.matrix(sampleMetadata)
  SampleSheet <- as.matrix(SampleSheet)
  write.table(sampleMetadata,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=T,quote=F,append=T,sep="\t")
  BamMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bam"]),"SampleName"],"Bam",sep="_"),SampleSheet[,"SampleName"])
  BigWigMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"bigwig"]),"SampleName"],"Bigwig",sep="_"),SampleSheet[,"SampleName"])
  IntervalMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[,"interval"]),"SampleName"],"Interval",sep="_"),SampleSheet[,"SampleName"])
  write.table("\n#sampleMapping",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("#Bams",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(BamMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("\n#BigWigs",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(BigWigMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table("\n#Intervals",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
  write.table(IntervalMappings,file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,append=T,sep="\t")
}
makeTrackHub <- function(QCobject,peaksDir,bigwigDir,IGVdirectory,genome){
  dir.create(IGVdirectory,showWarnings=F)  
  ss <- QCmetadata(QCobject)
  SampleSheet <- ss[,c("ID","Tissue","Factor")]
  colnames(SampleSheet)[1] <- "SampleName"
  fileSheet <- cbind(ss[,c("ID"),drop=F],rep(NA,nrow(ss)))
  peaks <- gsub("/+","/",dir(peaksDir,pattern="*peaks.bed",full.names=T))
  names(peaks) <- gsub("_WithInput.*","",basename(peaks))
  bigwigs <- gsub("/+","/",dir(bigwigDir,pattern="*.Dup.*\\.bw",full.names=T))
  names(bigwigs) <- gsub("DupMarkedNormalised.bw","",basename(bigwigs))
  fileSheet <- merge(fileSheet,cbind(names(peaks),peaks),all=T,by=1,all.x=T,all.y=F)
  fileSheet <- merge(fileSheet,cbind(names(bigwigs),bigwigs),all=T,by=1,all.x=T,all.y=F)
  colnames(SampleSheet) <- c("SampleName",'Tissue',"Factor")
  colnames(fileSheet) <- c("SampleName","bam","bigwig","interval") 
  MakeIGVSampleMetadata(SampleSheet,fileSheet,IGVdirectory)
  return(MakeIGVSessionXML(fileSheet,IGVdirectory,"IGVfull",genome,locusName="All"))
}

MakeIGVSessionXML <- function(SampleSheet,igvdirectory,XMLname,genomeName,locusName="All"){
  i <- 1
  require(XML)
  SampleSheet <- as.matrix(SampleSheet)
  Output <- file.path(igvdirectory,paste(XMLname,".xml",sep=""))
  GlobalNode <- newXMLNode("Global",attrs=c(genome.value=genomeName,groupTracksBy="Linking_id",locus=locusName,version=3))
  ResourcesNode <- newXMLNode("Resources",parent=GlobalNode)
  MetaDataNode <- newXMLNode("Resource",parent=ResourcesNode,attrs=c(name="SampleMetadata",path=relativePath(file.path(igvdirectory,"SampleMetadata.txt"),Output),relativePath=T))
  PanelDataNode <-  newXMLNode("Panel",attrs=c(height="350",name="DataPanel",width="1115"),parent=GlobalNode)
  #bamFiles <- SampleSheet[!is.na(SampleSheet[,"bam"]),"bam"]
  #bigwigFiles <- SampleSheet[!is.na(SampleSheet[,"bigwig"]),"bigwig"]
  #intervals <- SampleSheet[!is.na(SampleSheet[,"interval"]),"interval"]
  bamFiles <- SampleSheet[,"bam"]
  bigwigFiles <- SampleSheet[,"bigwig"]
  intervalFiles <- SampleSheet[,"interval"]    
  resources <- vector("list")
  #print(Output)
  for(i in 1:nrow(SampleSheet)){
    print(i)
    if(!is.na(SampleSheet[i,"bam"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Bam",sep="")
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bamFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",colorOption="UNEXPECTED_PAIR",displayMode="EXPANDED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bamFiles[i],Output),name=NewName,showDataRange="true",sortByTag="",visible="true"),parent=PanelDataNode)
    }
    if(!is.na(SampleSheet[i,"interval"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Interval",sep="")
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(intervalFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="45",id=relativePath(intervalFiles[i],Output),name=NewName,renderer="BASIC_FEATURE",showDataRange="true",sortable="false",visible="true",windowFunction="count"),parent=PanelDataNode)
    }
    if(!is.na(SampleSheet[i,"bigwig"])){
      NewName <- paste(SampleSheet[i,"SampleName"],"_Bigwig",sep="")
      print(relativePath(bigwigFiles[i],Output))
      resources <-  c(resources,list(newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName,name=NewName,path=relativePath(bigwigFiles[i],Output),relativePath=T))))
      TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoscale="true",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(bigwigFiles[i],Output),name=NewName,renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
      DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum="50",minimum="5",type="LINEAR"),parent=TrackNode)
    }
  }  
  saveXML(GlobalNode,file=Output)
  
  return(Output)
}

exportNormalisedBW <- function(bamFile,qc,normaliseTo="blacklisted"){
  require(GenomicAlignments)
  require(rtracklayer)
  #library(QuasR)
  extendBy <- fragmentlength(qc)
  #extendBy <- qc
  message("Reading tags from ",bamFile,appendLF=FALSE)
  #totalReads <- alignmentStats(bamFile)[,"mapped"]
  if(normaliseTo == "blacklisted"){
    totalReads <- qc@FlagAndTagCounts["MapQPass"] - qc@CountsInFeatures$BlackList
  }
  if(normaliseTo == "Total"){
    totalReads <- qc@FlagAndTagCounts["Mapped"]
  }  
  if(normaliseTo == "UniqueTotal"){
    totalReads <- qc@FlagAndTagCounts["Mapped"]-qc@FlagAndTagCounts["Duplicates"]
  }  
  total <- readGAlignmentsFromBam(bamFile)
  message("..done")
  message("Read in ",length(total)," reads")
  message("Extending reads to fragmentlength of ",extendBy," ..",appendLF=FALSE)
  temp <- resize(as(total,"GRanges"),extendBy,"start")
  message("..done")
  rm(total)
  gc()
  message("Calculating coverage..",appendLF=FALSE)
  genomeCov <- coverage(temp)
  rm(temp)
  message("..done")
  
  message("Normalised coverage..",appendLF=FALSE)
  genomeCov <- (genomeCov/totalReads)*1000000
  message("..done")
  message("Exporting coverage..",appendLF=FALSE)
  export.bw(genomeCov,file.path(dirname(dirname(bamFile)),gsub("\\.bam","Normalised\\.bw",basename(bamFile))))
  message("..done")
}

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
regionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,verbose=T,format="bam",seqlengths=NULL,forceFragment=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  result <- runRegionPlot(bamFile,testRanges,nOfWindows,FragmentLength,style,distanceAround,distanceInRegionStart,distanceOutRegionStart,distanceInRegionEnd,distanceOutRegionEnd,paired,normalize,plotBy,removeDup,format,seqlengths,forceFragment)
  return(result)  
}
runRegionPlot <- function(bamFile,testRanges,nOfWindows=100,FragmentLength=150,style="point",distanceAround=1500,distanceInRegionStart=1500,distanceOutRegionStart=1500,distanceInRegionEnd=1500,distanceOutRegionEnd=1500,paired=F,normalize="RPM",plotBy="coverage",removeDup=F,format="bam",seqlengths=NULL,forceFragment=NULL){
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
      #bamFile <- BamFile(bamFile,yieldSize=50000000L)
      #tempPaired <- readGAlignmentPairsFromBam(bamFile,param=Param)
      #tempPaired <- tempPaired[isProperPair(tempPaired)]
      gaPaired <- readGAlignmentsFromBam("/Users/tcarroll/Downloads//mergedETOH.bamRange5.bam", 
                                     param=ScanBamParam(what=c("mpos"),
                                     flag=scanBamFlag(isProperPair = TRUE,isFirstMateRead = TRUE)))      
      tempPos <- GRanges(seqnames(gaPaired[strand(gaPaired) == "+"]),
                      IRanges(
                        start=start(gaPaired[strand(gaPaired) == "+"]),
                        end=elementMetadata(gaPaired[strand(gaPaired) == "+"])$mpos
                        +qwidth(res3[strand(gaPaired) == "+"])))
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
    return(list(meansMat,AllRegionStart,AllRegionEnd,rbind(start,mid,end)))
    
  } 
}
