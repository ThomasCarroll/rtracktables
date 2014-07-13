MakeIGVSampleMetadata <- function(sampleMetadata,SampleSheet,igvdirectory){
  write.table("#sampleTable",file.path(igvdirectory,"SampleMetadata.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  colnames(sampleMetadata)[1] <- "Linking_id"
  print(colnames(sampleMetadata))
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
MakeIGVSessionXML <- function(files,fileType,names,igvdirectory,XMLname,genomeName,locusName="All"){
  i <- 1
  require(XML)
  files <- gsub(".xls",".bed",files)
  Output <- file.path(igvdirectory,paste(XMLname,".xml",sep=""))
  resources <- vector("list",length=length(files))
  if(fileType == "Bam"){
    NewName <- paste(names[i],"_Bam",sep="")
  }
  if(fileType == "Bigwig"){
    NewName <- paste(names[i],"_Bigwig",sep="")
  }
  if(fileType == "Macs"){
    NewName <- paste(names[i],"_Macs",sep="")
  }
  if(fileType == "MacsDiff"){
    NewName <- paste(names[i],"_MacsDiff",sep="")
  }
  if(fileType == "TPICs"){
    NewName <- paste(names[i],"_TPICs",sep="")
  }
  if(fileType == "Sicer"){
    NewName <- paste(names[i],"_Sicer",sep="")
  }
  GlobalNode <- newXMLNode("Global",attrs=c(genome.value=genomeName,groupTracksBy="Linking_id",locus=locusName,version=3))
  ResourcesNode <- newXMLNode("Resources",parent=GlobalNode)
  for(i in 1:length(resources)){
    resources[[i]] <-  newXMLNode("Resource",parent=ResourcesNode,attrs=c(label=NewName[i],name=NewName[i],path=relativePath(files[i],Output),relativePath=T))
  }
  MetaDataNode <- newXMLNode("Resource",parent=ResourcesNode,attrs=c(name="SampleMetadata",path=relativePath(file.path(igvdirectory,"SampleMetadata.txt"),Output),relativePath=T))
  PanelDataNode <-  newXMLNode("Panel",attrs=c(height="350",name="DataPanel",width="1115"),parent=GlobalNode)
  if(fileType == "Bam"){
    TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",colorOption="UNEXPECTED_PAIR",displayMode="EXPANDED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(files[i],Output),name=NewName[i],showDataRange="true",sortByTag="",visible="true"),parent=PanelDataNode)
  }
  if(fileType == "Bigwig"){
    TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",autoscale="true",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",id=relativePath(files[i],Output),name=NewName[i],renderer="BAR_CHART",showDataRange="true",visible="true",windowFunction="mean"),parent=PanelDataNode)
    DisplayRangeNode <-  newXMLNode("DataRange",attrs=c(baseline="0.0",drawBaseline="true",flipAxis="false",maximum="50",minimum="5",type="LINEAR"),parent=TrackNode)
  }
  if(fileType == "Macs" | fileType == "TPICs"| fileType == "Sicer" | fileType == "MacsDiff"){
    TrackNode <-  newXMLNode("Track",attrs=c(altColor="0,0,178",color="0,0,178",displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="45",id=relativePath(files[i],Output),name=NewName[i],renderer="BASIC_FEATURE",showDataRange="true",sortable="false",visible="true",windowFunction="count"),parent=PanelDataNode)
  }
  saveXML(GlobalNode,file=Output)
  
  return(Output)
}

SampleSheet <- cbind(
                  c("Pu1","Myc","Ik_prePro","Ik_pro"),
                  c("ProB","ProB","PreProB","ProB"),
                  c("Pu1","Myc","Ik","Ik")
)
colnames(SampleSheet) <- c("SampleName",'Tissue',"Factor")
fileSheet <- cbind(
  c("Pu1","Myc","Ik_prePro","Ik_pro"),
  c(NA,NA,NA,NA),
  c("Downloads/randomTracks-2//Pu1DupMarkedNormalised.bw",
    "Downloads/randomTracks-2//MycDupMarkedNormalised.bw",
    "Downloads/randomTracks-2//Ikaros_2_preproBDupMarkedNormalised.bw",
    "Downloads/randomTracks-2//Ikaros_1_proBDupMarkedNormalised.bw"),
  c("Downloads/randomTracks/Pu1_WithInput_Input_2_proB_peaks.bed",
    "Downloads/randomTracks/Myc_WithInput_Input_Ch12_peaks.bed",
    "Downloads/randomTracks/Ikaros_2_preproB_WithInput_Input_2_proB_peaks.bed",
    "Downloads/randomTracks/Ikaros_1_proB_WithInput_Input_2_proB_peaks.bed")
)
colnames(fileSheet) <- c("SampleName","bam","bigwig","interval") 
MakeIGVSampleMetadata(SampleSheet,fileSheet,getwd())
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
load(file ="Downloads/Sample_R2-INPUT6hDupMarked.RData")
bamFile <- "Downloads//Sample_R2-6hDupMarked.bam"
exportNormalisedBW(bamFile,qc)
?
gc
