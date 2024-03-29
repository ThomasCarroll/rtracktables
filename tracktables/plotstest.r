#' Plot regions
#'
#' A function to plot regions
#' 
#' @usage
#' \S4method{plotRegion}{ChIPprofile}(object)
#'
#'
#' @docType methods
#' @name plotRegion
#' @rdname plotRegion
#' @aliases plotRegion plotRegion,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @param object A ChIPprofile object 
#' @param gts A list of character vectors and GRanges
#' @param plotregion region to plot
plotRegion.ChIPprofile <- function(object,gts=NULL,plotregion="full")
{
  if(object@params$style=="region" & plotregion=="end"){
    bigTableFilt <- (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1))           
  }
  if(object@params$style=="region" & plotregion=="start"){
    bigTableFilt <- seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1))             
  } 
  #app <- lapply(gsets,function(x){colMeans(assays(object)[[1]][rowData(object)$name %in% x,])})
  if(!is.null(gts)){
    profileList <- list()
   for(p in 1:length(assays(object))){
     profileTemp <- assays(object)[[p]]
     if(object@params$style=="region" & (plotregion=="end" | plotregion=="start")){
     profileTempList <- lapply(gts,function(x)colMeans(profileTemp[rowData(object)$name %in% x,bigTableFilt])) 
     }else{
       profileTempList <- lapply(gts,function(x)colMeans(profileTemp[rowData(object)$name %in% x,]))        
     }
     profileMatTemp <- melt(as.data.frame(do.call(cbind,profileTempList)))
     if(object@params$style=="region" & plotregion=="full"){
     axisIndex=c(seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
                 (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                 (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                   seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)))
     }
     if(object@params$style=="point"){
       axisIndex=c(seq(1,(object@params$distanceAround+object@params$distanceAround+1)))
     }
     if(object@params$style=="region" & plotregion=="start"){
       axisIndex=bigTableFilt             
     }  
     if(object@params$style=="region" & plotregion=="end"){
       axisIndex=bigTableFilt             
     }       
     
     profileFrame <-data.frame("xIndex"=axisIndex,Group=profileMatTemp[,1],Sample=basename(unlist(exptData(object)["names"]))[p],Score=profileMatTemp[,2])

     profileList[[p]] <- profileFrame
   }  
   meltedProfileFrame <- do.call(rbind,profileList)
   colnames(meltedProfileFrame) <- c("xIndex","Group","Sample","Score")
  }else{
    profileList <- lapply(c(assays(object)),function(x)colMeans(x))
    if(object@params$style=="region" & (plotregion=="end" | plotregion=="start")){
      profileList <- lapply(c(assays(object)),function(x)colMeans(x[,bigTableFilt]))
    }else{
      profileFrame <- do.call(cbind,profileList)
    }
  }
    colnames(profileFrame) <- basename(unlist(exptData(object)["names"]))
    if(object@params$style=="region"){    
    axisIndex=c(seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                  seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)))
    }
    if(object@params$style=="point"){
      axisIndex=c(seq(1,(object@params$distanceAround+object@params$distanceAround+1)))
    }
    if(object@params$style=="region" & plotregion=="start"){
      axisIndex=bigTableFilt             
    }  
    if(object@params$style=="region" & plotregion=="end"){
      axisIndex=bigTableFilt             
    }       
    rownames(profileFrame) <- axisIndex
    meltedProfileFrame <- melt(profileFrame)
    colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")
  }
  #profileList <- lapply(c(assays(object)),function(y)lapply(gsets,function(x){colMeans(y[rowData(object)$name %in% x,])}))


  P <- ggplot(meltedProfileFrame,
              aes(x=xIndex,y=Score))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))
  if(object@params$style=="region" & plotregion=="full"){
  P <- P + scale_x_discrete(breaks=c(1,object@params$distanceInRegionStart+1,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1),
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+25*100,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+75*100,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100),
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+object@params$distanceInRegionStart+1,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
                       labels=c("Start-1500","Start","Start+1500","25%","75%","End-1500","End","End+1500"))+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
  if(object@params$style=="point"){
    P <- P + scale_x_discrete(breaks=c(1,object@params$distanceAround+1,object@params$distanceAround+1+object@params$distanceAround),
                              labels=c("Centre-1500","Centre","Centre+1500"))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
if(object@params$style=="region" & plotregion=="end"){
  P <- P + scale_x_discrete(breaks=seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)),
                            labels=c("End-1500","End","End+1500"))+
    theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
}
if(object@params$style=="region" & plotregion=="start"){
  P <- P + scale_x_discrete(breaks=seq(1,(object@params$distanceOutRegionstart+object@params$distanceInRegionStart+1)),
                            labels=c("End-1500","End","End+1500"))+
    theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
}  
  if(!is.null(gts)){
    P <- P+aes(group=Group)
  }
  return(P)
}

setGeneric("plotRegion", function(object="ChIPprofile",gts=NULL,plotregion="character") standardGeneric("plotRegion"))

#' @rdname plotRegion
#' @export
setMethod("plotRegion", signature(object="ChIPprofile"), plotRegion.ChIPprofile)

