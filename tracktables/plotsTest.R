#' Plot regions
#'
#' A function to plot regions
#' 
#' @usage
#' \S4method{plotRegion}{ChIPprofile}(object)
#'
#' @docType methods
#' @name plotRegion
#' @rdname plotRegion
#' @aliases plotRegion plotRegion,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @param object A ChIPprofile object 
plotRegion.ChIPprofile <- function(object,gsets=NULL)
{
  #app <- lapply(gsets,function(x){colMeans(assays(object)[[1]][rowData(object)$name %in% x,])})
  profileBy <- function(object,gsets){
    profileMat <- assays(object)[[1]]
    name = basename(unlist(exptData(object)["names"]))
    if(!is.null(gsets)){
    gsetsList <- lapply(gsets,function(x){colMeans(profileMat[rowData(object)$name %in% x,])})
    gsetsFrame <- do.call(cbind,gsetsList)
    colnames(gsetsFrame) <- names(gsets)
    }else{
      gsetsFrame <- cbind(colMeans(profileMat))
      colnames(gsetsFrame) <- "All"      
    }
    axisIndex=c(seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                  seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)))
    rownames(gsetsFrame) <- axisIndex
    meltedgsetsframe <- melt(gsetsFrame)
    meltedgsetsframe <- cbind(meltedgsetsframe,rep(name,nrow(meltedgsetsframe)))
    
  }
  for(i in 1:length(assays(object))){
    profileList[[i]] <-  profileBy(object,gsets)
  }
  #profileList <- lapply(c(assays(object)),function(y)lapply(gsets,function(x){colMeans(y[rowData(object)$name %in% x,])}))
  profileFrame <- do.call(function(x)rbind(x),profileList)
  colnames(profileFrame) <- c("xIndex","Group","Score","Sample")
  #profileFrame <- profileFrame[order(profileFrame[c("Sample","Group","xIndex")]),]
  P <- ggplot(profileFrame,
              aes(x=xIndex,y=Score,group=Sample,colour=Group))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))
  P <- P + scale_x_discrete(breaks=c(1,object@params$distanceInRegionStart+1,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1),
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+25*100,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+75*100,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100),
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+object@params$distanceInRegionStart+1,
                                (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
                       labels=c("Start-1500","Start","Start+1500","25%","75%","End-1500","End","End+1500"))+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  return(P)
}

setGeneric("plotRegion", function(object="ChIPprofile") standardGeneric("plotRegion"))

#' @rdname plotRegion
#' @export
setMethod("plotRegion", signature(object="ChIPprofile"), plotRegion.ChIPprofile)

