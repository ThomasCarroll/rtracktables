library(GenomicAlignments)
temp <- readGAlignmentsFromBam("/Users/tcarroll/Downloads/Sample_InpDN230614DupMarked.bam")
allchrs <- names(scanBamHeader("/Users/tcarroll/Downloads/Sample_InpDN230614DupMarked.bam")[[1]]$targets)
lengths <- as.vector(scanBamHeader("/Users/tcarroll/Downloads/Sample_InpDN230614DupMarked.bam")[[1]]$targets)
names(lengths) <- allchrs

tempCov <- coverage(temp)
tempCov3 <- tempCov[[10]]

TempBL <- reduce(as(slice(tempCov3,sd(tempCov3)*10),"IRanges"))
Distances <- end(TempBL[-1])-start(TempBL[-length(TempBL)])
plot(density(log10(Distances)),type="l",main="",xlab="log10_Min_Distance")
moreBL <- reduce(TempBL,min.gapwidth=1e+05)

save(tempCov,file="/home/pgellert/Dropbox
     (Lymphocyte_Developme)/tracktables/coverage.RData")


TempBL <- as(slice(tempCov3,sd(tempCov3)*10),"IRanges")
Distances <- end(TempBL[-1])-start(TempBL[-length(TempBL)])
plot(density(log10(Distances)),type="l")
plot(density(log2(viewSums(Views(tempCov3,TempBL))/width(TempBL))),xlim=c(3,5))

TempBL <- as(slice(tempCov3,median(sd(tempCov3))*1),"GRanges")

