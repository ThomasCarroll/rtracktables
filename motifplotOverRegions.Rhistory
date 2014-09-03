distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
distanceDownEndPos <- distanceUpStartPos
distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
distanceDownEndNeg <- distanceUpStartNeg
}else{
distanceUpStartPos <- distanceUpStart
distanceDownEndPos <- distanceDownEnd
distanceUpStartNeg <- distanceUpStart
distanceDownEndNeg <- distanceDownEnd
}
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
exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
)
message("...done")
style <- "percentOfRegion"
message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
testRangeNames <- unique(seqnames(testRanges))
temptestranges <- GRanges()
for(i in 1:length(testRangeNames)){
perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ])
#print(i)
perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ]
distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
}
distanceDownEnd <- distanceUpStart
}
message("..Done")
message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
testRanges <- temptestranges
message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
elementMetadata(testRanges) <- cbind(elementMetadata(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))
strand(testRanges[strand(testRanges) == "*"]) <- "+"
testRangesPos <- testRanges[strand(testRanges) == "+"]
testRangesNeg <- testRanges[strand(testRanges) == "-"]
if(style=="percentOfRegion"){
distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
distanceDownEndPos <- distanceUpStartPos
distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
distanceDownEndNeg <- distanceUpStartNeg
}else{
distanceUpStartPos <- distanceUpStart
distanceDownEndPos <- distanceDownEnd
distanceUpStartNeg <- distanceUpStart
distanceDownEndNeg <- distanceDownEnd
}
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
exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
)
message("...done")
if(style == "percentOfRegion"){
maxDistance <- round((distanceAround/100)*width(testRanges))
RegionsMat <- NULL
distanceUpStart <- NULL
distanceDownEnd <- NULL
}
maxDistance
if(is.null(seqlengths)){
seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
}else{
seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
}
lengths <- seqlengths(genomeCov)
allchrs <- names(lengths)
message("..Done")
style
if(style != "percentOfRegion"){
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
}
if(style == "percentOfRegion"){
message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
testRangeNames <- unique(seqnames(testRanges))
temptestranges <- GRanges()
for(i in 1:length(testRangeNames)){
perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ])
#print(i)
perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ]
distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
}
distanceDownEnd <- distanceUpStart
}
message("..Done")
message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
testRanges <- temptestranges
message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
elementMetadata(testRanges) <- cbind(elementMetadata(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))
strand(testRanges[strand(testRanges) == "*"]) <- "+"
testRangesPos <- testRanges[strand(testRanges) == "+"]
testRangesNeg <- testRanges[strand(testRanges) == "-"]
if(style=="percentOfRegion"){
distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
distanceDownEndPos <- distanceUpStartPos
distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
distanceDownEndNeg <- distanceUpStartNeg
}else{
distanceUpStartPos <- distanceUpStart
distanceDownEndPos <- distanceDownEnd
distanceUpStartNeg <- distanceUpStart
distanceDownEndNeg <- distanceDownEnd
}
message("..Done")
message("Found ",length(testRangesPos)," Watson strand regions")
message("Found ",length(testRangesNeg)," Crick strand regions")
## Extend regions and get positive versus negative reads.
message("Extending regions..",appendLF=FALSE)
exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
)
message("...done")
chromosomes <- seqlevels(genomeCov)
meansList <- vector("numeric")
grListWindowsPos <- GRanges()
grListWindowsNeg <- GRanges()
#grListWindows <- list()
message("Making windows..",appendLF=FALSE)
if(length(testRangesPos) > 0){
grWidths <- width(testRangesPos)+distanceUpStartPos+distanceDownEndPos
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
startPos <- start(testRangesPos)-((windows)*nOfWindows)-(windows/2)
for(i in 1:allWindows){
addToWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == allWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
(startPos)+(windows*(i-1)),
startPos+(windows*i)-1+addToWindow),giID=testRangesPos$giID))
}
grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
}
if(length(testRangesNeg) > 0){
grWidths <- width(testRangesNeg)
windows <- floor(grWidths%/%nOfWindows)
extraLastWindow <- grWidths%%nOfWindows
addToWindow <- 0
grWidths <- width(testRangesNeg)+distanceUpStartNeg+distanceDownEndNeg
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
endPos <- end(testRangesNeg)+(windows*nOfWindows)+round(windows/2)
for(i in 1:allWindows){
addToFirstWindow <- 0
addToLastWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == nOfWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
endPos-(windows*i)+1-addToWindow,
endPos-(windows*(i-1))),giID=testRangesNeg$giID))
}
grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
}
grListWindows <- c(grListWindowsPos,grListWindowsNeg)
message("..done\n")
for(c in 1:length(chromosomes)){
message("Processing inner region windows in ",chromosomes[c])
covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
doubleTemp <- viewMeans(covPerPeak)
names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$giID)
meansList <- c(meansList,doubleTemp)
message("..done")
message("Processing flanking windows in ",chromosomes[c])
tempstartRegionRangesPosMat <- NULL
tempendRegionRangesPosMat <- NULL
tempstartRegionRangesNegMat <- NULL
tempendRegionRangesNegMat <- NULL
}
meansMat <- matrix(meansList,ncol=allWindows,byrow=T)
rownames(meansMat) <- matrix(names(meansList),ncol=allWindows,byrow=T)[,1]
plot(apply(meansMat,2,function(x)mean(x,na.rm = T)),type="l")
motif75 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="75%",removeRand=T)
motif80 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="80%",removeRand=T)
motif85 <- pwmToCoverage(mdb.ctcf[[1]],Mmusculus,min="85%",removeRand=T)
maxDistance <- round((distanceAround/100)*width(testRanges))
RegionsMat <- NULL
distanceUpStart <- NULL
distanceDownEnd <- NULL
genomeCov <- motif75
if(is.null(seqlengths)){
seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
}else{
seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
}
lengths <- seqlengths(genomeCov)
allchrs <- names(lengths)
message("..Done")
testRanges <- DP_thymocytes
if(style == "percentOfRegion"){
message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
testRangeNames <- unique(seqnames(testRanges))
temptestranges <- GRanges()
for(i in 1:length(testRangeNames)){
perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ])
#print(i)
perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ]
distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
}
distanceDownEnd <- distanceUpStart
}
message("..Done")
message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
testRanges <- temptestranges
message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
elementMetadata(testRanges) <- cbind(elementMetadata(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))
strand(testRanges[strand(testRanges) == "*"]) <- "+"
testRangesPos <- testRanges[strand(testRanges) == "+"]
testRangesNeg <- testRanges[strand(testRanges) == "-"]
if(style=="percentOfRegion"){
distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
distanceDownEndPos <- distanceUpStartPos
distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
distanceDownEndNeg <- distanceUpStartNeg
}else{
distanceUpStartPos <- distanceUpStart
distanceDownEndPos <- distanceDownEnd
distanceUpStartNeg <- distanceUpStart
distanceDownEndNeg <- distanceDownEnd
}
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
exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
)
message("...done")
testRanges
testRanges <- DP_thymocytes
testRangeNames
testRanges
maxDistance <- round((distanceAround/100)*width(testRanges))
RegionsMat <- NULL
distanceUpStart <- NULL
distanceDownEnd <- NULL
genomeCov <- bamFile
if(is.null(seqlengths)){
seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
}else{
seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
}
lengths <- seqlengths(genomeCov)
allchrs <- names(lengths)
message("..Done")
genomeCov <- motif75
if(is.null(seqlengths)){
seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
}else{
seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
}
lengths <- seqlengths(genomeCov)
allchrs <- names(lengths)
message("..Done")
style
if(style == "percentOfRegion"){
message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
testRangeNames <- unique(seqnames(testRanges))
temptestranges <- GRanges()
for(i in 1:length(testRangeNames)){
perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ])
#print(i)
perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
& start(perchrRanges)-perChrMaxDistance > 0 ]
distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
}
distanceDownEnd <- distanceUpStart
}
message("..Done")
message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
testRanges <- temptestranges
message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
elementMetadata(testRanges) <- cbind(elementMetadata(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))
strand(testRanges[strand(testRanges) == "*"]) <- "+"
testRangesPos <- testRanges[strand(testRanges) == "+"]
testRangesNeg <- testRanges[strand(testRanges) == "-"]
if(style=="percentOfRegion"){
distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")]
distanceDownEndPos <- distanceUpStartPos
distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
distanceDownEndNeg <- distanceUpStartNeg
}else{
distanceUpStartPos <- distanceUpStart
distanceDownEndPos <- distanceDownEnd
distanceUpStartNeg <- distanceUpStart
distanceDownEndNeg <- distanceDownEnd
}
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
exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
)
message("...done")
meansList <- vector("numeric")
grListWindowsPos <- GRanges()
grListWindowsNeg <- GRanges()
#grListWindows <- list()
message("Making windows..",appendLF=FALSE)
if(length(testRangesPos) > 0){
grWidths <- width(testRangesPos)+distanceUpStartPos+distanceDownEndPos
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
startPos <- start(testRangesPos)-((windows)*nOfWindows)-(windows/2)
for(i in 1:allWindows){
addToWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == allWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
(startPos)+(windows*(i-1)),
startPos+(windows*i)-1+addToWindow),giID=testRangesPos$giID))
}
grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
}
if(length(testRangesNeg) > 0){
grWidths <- width(testRangesNeg)
windows <- floor(grWidths%/%nOfWindows)
extraLastWindow <- grWidths%%nOfWindows
addToWindow <- 0
grWidths <- width(testRangesNeg)+distanceUpStartNeg+distanceDownEndNeg
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
endPos <- end(testRangesNeg)+(windows*nOfWindows)+round(windows/2)
for(i in 1:allWindows){
addToFirstWindow <- 0
addToLastWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == nOfWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
endPos-(windows*i)+1-addToWindow,
endPos-(windows*(i-1))),giID=testRangesNeg$giID))
}
grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
}
grListWindows <- c(grListWindowsPos,grListWindowsNeg)
message("..done\n")
for(c in 1:length(chromosomes)){
message("Processing inner region windows in ",chromosomes[c])
covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
doubleTemp <- viewMeans(covPerPeak)
names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$giID)
meansList <- c(meansList,doubleTemp)
message("..done")
message("Processing flanking windows in ",chromosomes[c])
tempstartRegionRangesPosMat <- NULL
tempendRegionRangesPosMat <- NULL
tempstartRegionRangesNegMat <- NULL
tempendRegionRangesNegMat <- NULL
}
meansMat <- matrix(meansList,ncol=allWindows,byrow=T)
rownames(meansMat) <- matrix(names(meansList),ncol=allWindows,byrow=T)[,1]
plot(apply(meansMat,2,function(x)mean(x,na.rm = T)),type="l")
genomeCov <- motif80
meansList <- vector("numeric")
grListWindowsPos <- GRanges()
grListWindowsNeg <- GRanges()
#grListWindows <- list()
message("Making windows..",appendLF=FALSE)
if(length(testRangesPos) > 0){
grWidths <- width(testRangesPos)+distanceUpStartPos+distanceDownEndPos
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
startPos <- start(testRangesPos)-((windows)*nOfWindows)-(windows/2)
for(i in 1:allWindows){
addToWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == allWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(
(startPos)+(windows*(i-1)),
startPos+(windows*i)-1+addToWindow),giID=testRangesPos$giID))
}
grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
}
if(length(testRangesNeg) > 0){
grWidths <- width(testRangesNeg)
windows <- floor(grWidths%/%nOfWindows)
extraLastWindow <- grWidths%%nOfWindows
addToWindow <- 0
grWidths <- width(testRangesNeg)+distanceUpStartNeg+distanceDownEndNeg
allWindows <- nOfWindows*3
windows <- floor(grWidths%/%allWindows)
extraLastWindow <- grWidths%%allWindows
addToWindow <- 0
endPos <- end(testRangesNeg)+(windows*nOfWindows)+round(windows/2)
for(i in 1:allWindows){
addToFirstWindow <- 0
addToLastWindow <- 0
if(i == 1){
addToWindow <- round(extraLastWindow/2)
}
if(i == nOfWindows){
addToWindow <- round(extraLastWindow/2)
}
grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(
endPos-(windows*i)+1-addToWindow,
endPos-(windows*(i-1))),giID=testRangesNeg$giID))
}
grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
}
grListWindows <- c(grListWindowsPos,grListWindowsNeg)
message("..done\n")
for(c in 1:length(chromosomes)){
message("Processing inner region windows in ",chromosomes[c])
covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
doubleTemp <- viewMeans(covPerPeak)
names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$giID)
meansList <- c(meansList,doubleTemp)
message("..done")
message("Processing flanking windows in ",chromosomes[c])
tempstartRegionRangesPosMat <- NULL
tempendRegionRangesPosMat <- NULL
tempstartRegionRangesNegMat <- NULL
tempendRegionRangesNegMat <- NULL
}
meansMat <- matrix(meansList,ncol=allWindows,byrow=T)
rownames(meansMat) <- matrix(names(meansList),ncol=allWindows,byrow=T)[,1]
plot(apply(meansMat,2,function(x)mean(x,na.rm = T)),type="l")
plot(caTools::runmean(apply(meansMat,2,function(x)mean(x,na.rm = T)),50),type="l")
plot(caTools::runmean(apply(meansMat,2,function(x)mean(x,na.rm = T)),10),type="l")
plot(caTools::runmean(apply(meansMat,2,function(x)mean(x,na.rm = T)),30),type="l")
plot(caTools::runmean(apply(meansMat,2,function(x)mean(x,na.rm = T)),20),type="l")
savehistory("~/Documents/t2/rtracktables/motifplotOverRegions.Rhistory")
