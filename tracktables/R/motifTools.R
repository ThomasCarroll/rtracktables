#' PWM hits as an RLElist
#'
#' Creates rlelist of pwm hits. 
#' 
#'
#'
#' @docType methods
#' @name pwmToCoverage
#' @rdname pwmToCoverage
#' 
#' @author Thomas Carroll
#'
#' @param pwm A pwn matrix object.
#' @param genome A BSgenome object
#' @param min pwm score (as percentage of maximum score) cutoff
#' @param removeRand Remove contigs with rand string
#' @export
pwmToCoverage <- function(pwm,genome,min="70%",removeRand=FALSE){
  
  allchrs <- seqnames(genome)
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  intergerList <- lapply(allchrs,function(x)pwmHitAsCoverage(pwm,genome,min,x))
  myrle <- RleList(intergerList,compress=F)
  names(myrle) <- allchrs
  myrle
}  

pwmHitAsCoverage <- function(pwm,genome,min,chrofinterest){
  posMotifs <- matchPWM(pwm,genome[[chrofinterest]],min.score=min)
  negMotifs <- matchPWM(reverseComplement(pwm),genome[[chrofinterest]],min.score=min)
  if(length(posMotifs) > 0){
    rleMotifHitPos <- coverage(GRanges(seqnames=chrofinterest,ranges(posMotifs),strand="+"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitPos <- RleList(rep(0,length(genome[[chrofinterest]])))
  }
  if(length(negMotifs) > 0){
    rleMotifHitNeg <- coverage(GRanges(seqnames=chrofinterest,ranges(negMotifs),strand="-"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitNeg <- RleList(rep(0,length(genome[[chrofinterest]])))    
  }
  rleTotal <- rleMotifHitPos+rleMotifHitNeg
  return(rleTotal[[1]])
}

pwmHitAsGRanges <- function(pwm,genome,min,chrofinterest){
  posMotifs <- matchPWM(pwm,genome[[chrofinterest]],min.score=min,with.score=T)
  negMotifs <- matchPWM(reverseComplement(pwm),genome[[chrofinterest]],min.score=min,with.score=T)
  posMotifs <- GRanges(seqnames=rep(chrofinterest,length(posMotifs)),
                       ranges=ranges(posMotifs),
                       strand=rep("+",length(posMotifs))                      
                       )
  negMotifs <- GRanges(seqnames=rep(chrofinterest,length(negMotifs)),
                       ranges=ranges(negMotifs),
                       strand=rep("-",length(negMotifs))                      
                        )
  #strand(negMotifs) <- rep("-",length(negMotifs))
  GRangesTotal <- c(posMotifs,negMotifs)
  return(GRangesTotal)
}

pwmToGranges <- function(pwm,genome,min,chrs=NULL){
  if(is.null(chrs)){
    chrs <- seqnames(Mmusculus)
  }
  chrs <- unique(chrs)
  Res <- lapply(chrs,function(x)pwmHitAsGRanges(pwm,genome,min,x))
  pwmHitsAsGRanges <- unlist(GRangesList(unlist(Res)))
}

makeGRangesWithSummary <- function(GRangesSet,scoreBy){
  temp <- viewSums(Views(scoreBy,
                 ranges(GRangesSet)
            )
  )  
  elementMetadata(GRangesSet) <- data.frame(score=temp)
  return(GRangesSet)
}
  

rleFromScoresInGRanges <- function(GRangesSet,scoreBy,chrs){
  chrs <- chrs[chrs %in% names(scoreBy)]
  chrs <- chrs[chrs %in% names(seqlengths(GRangesSet))]  
  res <- lapply(chrs,function(x)
                makeGRangesWithSummary(
                  GRangesSet[
                    seqnames(GRangesSet) %in% x,],
                  scoreBy[[x]]        
                  ))
  return(unlist(GRangesList(unlist(res))))
}

motifCov <- function(genome,regions,pwm,chrOfInterest,atCentre=FALSE){
  reducedregions <- reduce(regions[seqnames(regions) %in% chrOfInterest])
  regionViews <- Views(genome[[chrOfInterest]],ranges(reducedregions))
  trial <- matchPWM(pwm,regionViews,min.score = 0,with.score = T)
  if(atCentre==TRUE){
    theRanges <- resize(as(trial,"IRanges"),1,"center")
  }
  if(atCentre==FALSE){
    theRanges <- as(trial,"IRanges")
  }
  if(length(theRanges) > 0){
    motifCov <- unlist(coverage(GRanges(chrOfInterest,theRanges,"*",elementMetadata(trial)$score),weight="elementMetadata.trial..score"))
  }else{
    motifCov <- unlist(RleList(rep(0,length(genome[[chrOfInterest]]))))
  }
  return(motifCov)
}

#' Motif score as an RLElist
#'
#' Creates rlelist of pwm scores 
#'
#'
#' @docType methods
#' @name makeMotifScoreRle
#' @rdname makeMotifScoreRle
#' 
#' @author Thomas Carroll
#'
#' @param pwm A pwn matrix object.
#' @param regions GRanges object to include in pwm rlelist
#' @param genome A BSgenome object
#' @param extend bps to extend regions by
#' @param removeRand Remove contigs with rand string
#' @param strandScore Method for averaging strand. Options are max, mean, sum, bothstrands
#' @export
makeMotifScoreRle <- function(pwm,regions,genome,extend,removeRand=FALSE,strandScore="mean",atCentre=FALSE){
  regions <- GRanges(seqnames(regions),IRanges(start(regions)-extend,end(regions)+extend),strand=strand(regions),elementMetadata(regions))
  lengths <- seqlengths(genome)
  ## Filter testRanges to those contained within chromosomes.
  message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
  testRangeNames <- unique(seqnames(regions))
  temptestranges <- GRanges()
  maxDistance <- extend
  for(i in 1:length(testRangeNames)){
    perchrRanges <- regions[seqnames(regions) %in% as.vector(testRangeNames[i])]
    temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                    & start(perchrRanges)-maxDistance > 0 ])
    #print(i)
  }
  
  message("..Done")
  message("Filtered ",length(regions)-length(temptestranges)," of ",length(regions)," regions")
  
  regions <- temptestranges
  
  allchrs <- as.vector(unique(seqnames(regions)))
  
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  forMotif <- list()
  revMotif <- list()
  message("Scoring motifs on positive strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])
    forMotif[[k]] <- unlist(motifCov(genome,regions,pwm,allchrs[k],atCentre))  
  }
  message("..done")
  
  message("Scoring motifs on negative strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])    
    revMotif[[k]] <- unlist(motifCov(genome,regions,reverseComplement(pwm),allchrs[k]))    
  }
  message("..done")
  
  #motifScoreRLEComplement <- lapply(allchrs,function(x)motifCov(genome,regions,reverseComplement(pwm),x))
  #myrle <- RleList(motifScoreRLE,compress=F)
  # myrleComplement <- RleList(motifScoreRLEComplement,compress=F)
  revMotif <- RleList(revMotif,compress=F)
  forMotif <- RleList(forMotif,compress=F)
  if(strandScore=="sum"){
    MotifScore <- revMotif+forMotif
    names(MotifScore) <- allchrs    
  }
  if(strandScore=="mean"){  
    MotifScore <- (revMotif+forMotif)/2
    names(MotifScore) <- allchrs    
  }
  if(strandScore=="max"){  
    MotifScore <- pmax(revMotif,forMotif)
    names(MotifScore) <- allchrs
  }  
  if(strandScore=="bothstrands"){  
    MotifScore <- list(revMotif,forMotif)
    names(MotifScore[[1]]) <- allchrs
    names(MotifScore[[2]]) <- allchrs    
  }  
  return(MotifScore)
}
