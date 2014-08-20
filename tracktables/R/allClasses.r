#' The ChIPprofile class
#' @export
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(params="list"
         ))

#' Merge ChIPprofile objects
#'
#' A Merge two ChIPprofiles into 1 combined ChIPprofile
#' 
#' @usage
#' \S4method{mergeChIPprofiles}{ChIPprofile,ChIPprofile}(ChIPprofile)
#'
#' @docType methods
#' @name mergeChIPprofiles
#' @rdname mergeChIPprofiles
#' @aliases mergeChIPprofiles mergeChIPprofiles,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @export
#' @param x A ChIPprofile object 
#' @param y A ChIPprofile object
mergeChIPprofiles <- function(x,y){
  t1 <- rowData(x)
  filtDupt1 <- as.numeric(names(table(as.matrix(findOverlaps(t1,t1,type="equal"))[,2])[
    table(as.matrix(findOverlaps(t1,t1,type="equal"))[,2]) >1]))
  if(length(filtDupt1)){
    t1assays <- lapply(assays(x),function(x)x[-filtDupt1,])
    t1 <- t1[-filtDupt1,]
  }else{
    t1assays <- assays(x)
  }
  t2 <- rowData(y)
  filtDupt2 <- as.numeric(names(table(as.matrix(findOverlaps(t2,t2,type="equal"))[,2])[
    table(as.matrix(findOverlaps(t2,t2,type="equal"))[,2]) >1]))
  if(length(filtDupt2)){
    t2assays <- lapply(assays(y),function(x)x[-filtDupt2,])
    t2 <- t2[-filtDupt2,]
  }else{
    t1assays <- assays(x)
  }
  t3 <- t1[as.matrix(findOverlaps(t1,t2,type="equal"))[,2]]
  t4 <- t2[as.matrix(findOverlaps(t1,t2,type="equal"))[,1]]
  
  bothMeta <- NULL
  if(!is.null(elementMetadata(t3)) & !is.null(elementMetadata(t4))){
    tempId <- seq(1,length(findOverlaps(t1,t2,type="equal"))) 
    t3Data <- cbind(tempId,as.data.frame(elementMetadata(t3)))
    t4Data <- as.data.frame(elementMetadata(t4))
    # Can fix this to always get some metadata if any grange does
    t4Data <- cbind(tempId,t4Data[,!colnames(t4Data) %in% colnames(t3Data)])
    bothMeta <- merge(t3Data,t4Data,by="tempId",all=T)    
  }
  elementMetadata(t3) <- bothMeta[,!colnames(bothMeta) %in% "tempId"]
  filtIndexX <- as.matrix(findOverlaps(t1,t2,type="equal"))[,2]
  filtIndexY <- as.matrix(findOverlaps(t1,t2,type="equal"))[,1]
  filtAssayListX <- lapply(t1assays,function(x)x[filtIndexX,])
  filtAssayListY <- lapply(t2assays,function(x)x[filtIndexY,])
  profileSample <- SummarizedExperiment(c(filtAssayListX,filtAssayListY),rowData=t3)
  exptData(profileSample) <- list(names=c(exptData(x)$names,exptData(y)$names))
  paramList <- x@params
  return(new("ChIPprofile",profileSample,params=paramList))
  
}

setGeneric("mergeChIPprofiles", function(x="ChIPprofile",y="ChIPprofile") standardGeneric("mergeChIPprofiles"))

