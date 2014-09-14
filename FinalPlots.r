library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)
library(MotifDb)
library(Biostrings)
library(seqLogo)
#load("/home/pgellert/MatthiasTrial/superenhancers.RData")
load("/Users/tcarroll/Downloads/superenhancers.RData")

#dpctcfPeaks <- ChIPQC:::GetGRanges("/home//pgellert/Dropbox (Lymphocyte_Developme)/tracktables/DP_CTCF_WithInput_DP_Input_peaks.bed")
#dpctcfPeaks <- ChIPQC:::GetGRanges("/Users/tcarroll/Downloads/DP_CTCF_WithInput_DP_Input_peaks.bed")


mdb.ctcf <- MotifDb [grep ('ctcf', values (MotifDb)$geneSymbol, ignore.case=TRUE)]


files <- dir("/Users/tcarroll/Downloads/VladSeitan/",pattern="*.bed",full.names=T)
SElist <- lapply(files,function(x)ChIPQC:::GetGRanges(x,simplify=T))
SElist <- c(list(DP_thymocytes),list(ESCs),SElist)
names(SElist) <- c("Ing_DP","Ing_ESC",gsub("\\.bed\\.csv","",basename(files)))

SElist2 <- lapply(SElist,function(x)GRanges(seqnames(x),IRanges(start(x)-(width(x)*2),end(x)+(width(x)*2)),strand="+"))

extendAllSE <- reduce(unlist(GRangesList(SElist2)))
motifScores_DP_thy_Enh_Max4 <- makeMotifScoreRle(mdb.ctcf[[1]],extendAllSE,Mmusculus,1000,removeRand=TRUE,strandScore="max")
SElist[[1]] <- ChIPQC:::GetGRanges(SElist[[1]],simplify=T)
SElist[[2]] <- ChIPQC:::GetGRanges(SElist[[2]],simplify=T)
allSE <- unlist(GRangesList(SElist))
elementMetadata(allSE) <- data.frame(name=names(allSE))

ctcf_all_motifNew <- regionPlot(motifScores_DP_thy_Enh_Max4,allSE,nOfWindows=100,style="percentOfRegion",format="rlelist",FragmentLength=130,distanceAround = 100)


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))[-c(1,2)]))

names(gts) <- c(unique(names(allSE)),"All_SE")
library(ggplot2)
p <- plotRegion(ctcf_all_motifNew,gts=gts)
p <- p+facet_wrap(~Group)+aes(colour=Group)







motifScores_DP_thy_Enh_Max4 <- makeMotifScoreRle(mdb.ctcf[[1]],extendAllSE,Mmusculus,1000,removeRand=TRUE,strandScore="max")
SElist[[1]] <- ChIPQC:::GetGRanges(SElist[[1]],simplify=T)
SElist[[2]] <- ChIPQC:::GetGRanges(SElist[[2]],simplify=T)
allSE <- unlist(GRangesList(SElist))
elementMetadata(allSE) <- data.frame(name=names(allSE))



ctcf_all_motifNew <- regionPlot(motifScores_DP_thy_Enh_Max4,allSE,nOfWindows=100,style="percentOfRegion",format="rlelist",FragmentLength=130,distanceAround = 100)


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))[-c(1,2)]))

names(gts) <- c(unique(names(allSE)),"All_SE")

p <- plotRegion(ctcf_all_motifNew,gts=gts)
p <- p+facet_wrap(~Group)+aes(colour=Group)

motifScores_DP_thy_Enh_Max4 <- makeMotifScoreRle(mdb.ctcf[[1]],extendAllSE,Mmusculus,1000,removeRand=TRUE,strandScore="max")
SElist[[1]] <- ChIPQC:::GetGRanges(SElist[[1]],simplify=T)
SElist[[2]] <- ChIPQC:::GetGRanges(SElist[[2]],simplify=T)
allSE <- unlist(GRangesList(SElist))
elementMetadata(allSE) <- data.frame(name=names(allSE))

ctcf_all_motifNew <- regionPlot(motifScores_DP_thy_Enh_Max4,allSE,nOfWindows=100,style="percentOfRegion",format="rlelist",FragmentLength=130,distanceAround = 100)


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))[-c(1,2)]))

names(gts) <- c(unique(names(allSE)),"All_SE")
library(ggplot2)
p <- plotRegion(ctcf_all_motifNew,gts=gts)
p <- p+facet_wrap(~Group)+aes(colour=Group)



SElist2 <- lapply(SElist,function(x)GRanges(seqnames(x),IRanges(start(x)-(width(x)*2),end(x)+(width(x)*2)),strand="+"))

extendAllSE <- reduce(unlist(GRangesList(SElist2[c(1,7)])))
motifScores_DP_thy_Enh_Max4 <- makeMotifScoreRle(mdb.ctcf[[1]],extendAllSE,Mmusculus,1000,removeRand=TRUE,strandScore="max",atCentre = T)

SElist[[1]] <- ChIPQC:::GetGRanges(SElist[[1]],simplify=T)
SElist[[2]] <- ChIPQC:::GetGRanges(SElist[[2]],simplify=T)
allSE <- unlist(GRangesList(SElist[c(1,7)]))
elementMetadata(allSE) <- data.frame(name=names(allSE))

ctcf_all_motifNew <- regionPlot(motifScores_DP_thy_Enh_Max4,allSE,nOfWindows=100,style="percentOfRegion",format="rlelist",FragmentLength=130,distanceAround = 100)


gts <- c(as.list(unique(names(allSE))),list(unique(names(allSE))))

names(gts) <- c(unique(names(allSE)),"All_SE")

p <- plotRegion(ctcf_all_motifNew,gts=gts)
p <- p+facet_wrap(~Group)+aes(colour=Group)

temp <- matrix(as.integer(mdb.ctcf[[1]]*10^9),nrow=4,byrow=F,ncol=19)
rownames(temp) <- c("A","C","G","T")

genomic.acgt = getBackgroundFrequencies(Mmusculus)
pwms.denovo = PFMtoPWM(temp, prior=genomic.acgt)

ctcfVsFlank <- function(x){
  gt <- x
  gt <- gt[width(gt) > 1021]
  
  Edges<- c(
    GRanges(seqnames(gt),IRanges(start(gt)-500,start(gt)+500),strand="+"),
    GRanges(seqnames(gt),IRanges(end(gt)-500,end(gt)+500),strand="+"))
  
  fl <- c(
    GRanges(seqnames(gt),IRanges(start(gt)-500-width(gt),start(gt)-500),strand="+"),
    GRanges(seqnames(gt),IRanges(end(gt)+500,end(gt)+500+width(gt)),strand="+"))
  
  EdgesSeq <- getSeq(Mmusculus,Edges)
  flSeq <- getSeq(Mmusculus,fl)
  tempBack <- makePWMLognBackground(flSeq,pwms.denovo)  
  res.EdgesSeq = motifEnrichment(EdgesSeq,tempBack,group.only=T)
  EdgesSeqRes <- groupReport(res.EdgesSeq)
  TIngRes <- EdgesSeqRes
  TIngRes
}
resAllvsFlank <- lapply(SElist[1],ctcfVsFlank)


