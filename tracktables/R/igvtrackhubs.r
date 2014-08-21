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
