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
plotRegion.ChIPprofile <- function(object)
{
  profileList <- lapply(c(assays(object)),colMeans)
  profileFrame <- do.call(cbind,profileList)
  colnames(profileFrame) <- basename(unlist(exptData(object)["names"]))
  axisIndex=c(seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)),
              (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
              (object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                seq(1,(object@params$distanceInRegionStart+object@params$distanceInRegionStart+1)))
  rownames(profileFrame) <- axisIndex
  meltedProfileFrame <- melt(profileFrame)
  colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")
  P <- ggplot(meltedProfileFrame,
              aes(x=xIndex,y=Score))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))
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

