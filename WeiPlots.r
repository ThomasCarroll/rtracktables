
preB1 <- read.delim("/Users/tcarroll//Downloads/GSM1296547_12078.rpkm.txt",sep="\t",h=T)
preB2 <- read.delim("/Users/tcarroll//Downloads/GSM1296548_12079.rpkm.txt",sep="\t",h=T)
preB1 <- cbind(gsub(":.*","",preB1[,1]),preB1)
preB2 <- cbind(gsub(":.*","",preB2[,1]),preB2)

mergedPreB <- merge(preB1,preB2,by=1,all=T)
mergedPreBRes <- cbind(as.vector(mergedPreB[,1]),apply(mergedPreB,1,function(x) as.numeric(x[3])+as.numeric(x[5])/2))

mergedPreBRes <- mergedPreBRes[mergedPreBRes[,1] %in% toupper(mm9PC$name),]
mergedPreBRes <- mergedPreBRes[order(as.numeric(mergedPreBRes[,2]),decreasing=T),]
top1000byRNAseq <- capitalize(tolower(mergedPreBRes[1:1000,1]))
mergedPreBRes <- mergedPreBRes[order(as.numeric(mergedPreBRes[,2]),decreasing=F),]
bottom1000byRNAseq <- capitalize(tolower(mergedPreBRes[1:1000,1]))

nuc4Top <- c(colMeans(nuc0hrange4[[1]][rownames(nuc0hrange4[[3]]) %in% top1000byRNAseq,]),
            colMeans(nuc0hrange4[[2]][rownames(nuc0hrange4[[3]]) %in% top1000byRNAseq,]),
            colMeans(nuc0hrange4[[3]][rownames(nuc0hrange4[[3]]) %in% top1000byRNAseq,])
)

nucTop <- c(colMeans(nuc0hrange5[[1]][rownames(nuc0hrange5[[3]]) %in% top1000byRNAseq,]),
colMeans(nuc0hrange5[[2]][rownames(nuc0hrange5[[3]]) %in% top1000byRNAseq,]),
colMeans(nuc0hrange5[[3]][rownames(nuc0hrange5[[3]]) %in% top1000byRNAseq,])
)
plot(nuc0hrange5[[4]][,1],nucTop,type="l",lwd=2)
plot(nuc0hrange5[[4]][1:3000,1],runmean(nuc0hrange5[[4]][1:3000,2],1),col="darkred",type="l")

polTop <- c(colMeans(pol0h[[1]][rownames(pol0h[[3]]) %in% top1000byRNAseq,]),
            colMeans(pol0h[[2]][rownames(pol0h[[3]]) %in% top1000byRNAseq,]),
            colMeans(pol0h[[3]][rownames(pol0h[[3]]) %in% top1000byRNAseq,])
)
plot(pol0h[[4]][,1],runmean(polTop,20),type="l",col="darkblue",lwd=2)
lines(pol0h[[4]][,1],runmean(polBottom,20),col="darkred",lwd=2)
lines(pol0h[[4]][,1],runmean(pol0h[[4]][,2],10),col="darkgreen")
plot(nuc0hrange5[[4]][1:3000,1],nucBottom[1:3000],type="l",lwd=2)
png("nuctop.png")
plot(nuc0hrange5[[4]][1:3000,1],nucTop[1:3000],col="darkred",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nucbottom.png")
plot(nuc0hrange5[[4]][1:3000,1],nucBottom[1:3000],col="darkblue",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nuc.png")
plot(nuc0hrange5[[4]][1:3000,1],nuc0hrange5[[4]][1:3000,2],lwd=1,ylim=c(0,5),type="l")
dev.off()

png("nucRange4top.png")
plot(nuc0hrange4[[4]][1:3000,1],nuc4Top[1:3000],col="darkred",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nucRange4bottom.png")
plot(nuc0hrange4[[4]][1:3000,1],nuc4Bottom[1:3000],col="darkblue",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nucRange4.png")
plot(nuc0hrange4[[4]][1:3000,1],nuc0hrange4[[4]][1:3000,2],lwd=1,ylim=c(0,5),type="l")
dev.off()

png("poltop.png")
plot(pol0h[[4]][1:3000,1],polTop[1:3000],col="darkred",lwd=1,ylim=c(0,10),type="l")
dev.off()
png("polbottom.png")
plot(pol0h[[4]][1:3000,1],polBottom[1:3000],col="darkblue",lwd=1,ylim=c(0,10),type="l")
dev.off()
png("pol.png")
plot(pol0h[[4]][1:3000,1],runmean(pol0h[[4]][1:3000,2],10),lwd=1,ylim=c(0,10),type="l")
dev.off()


png("nuctopWhole.png")
plot(nuc0hrange5[[4]][,1],nucTop,col="darkred",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nucbottomWhole.png")
plot(nuc0hrange5[[4]][,1],nucBottom,col="darkblue",lwd=1,ylim=c(0,5),type="l")
dev.off()
png("nucWhole.png")
plot(nuc0hrange5[[4]][,1],nuc0hrange5[[4]][,2],lwd=1,ylim=c(0,5),type="l")
dev.off()

png("poltopWhole.png")
plot(pol0h[[4]][,1],polTop,type="l",col="darkred",lwd=1,ylim=c(0,10))
dev.off()
png("polbottomWhole.png")
plot(pol0h[[4]][,1],polBottom,col="darkblue",lwd=1,ylim=c(0,10),type="l")
dev.off()
png("polWhole.png")
plot(pol0h[[4]][,1],runmean(pol0h[[4]][,2],10),lwd=1,ylim=c(0,10),type="l")
dev.off()

plot(pol0h[[4]][1:3000,1],runmean(pol0h[[4]][1:3000,2],10),col="darkred",type="l")

polBottom <- c(colMeans(pol0h[[1]][rownames(pol0h[[3]]) %in% bottom1000byRNAseq,]),
            colMeans(pol0h[[2]][rownames(pol0h[[3]]) %in% bottom1000byRNAseq,]),
            colMeans(pol0h[[3]][rownames(pol0h[[3]]) %in% bottom1000byRNAseq,])
)
plot(pol0h[[4]][,1],polBottom,type="l",lwd=2)
plot(pol0h[[4]][1:3000,1],runmean(pol0h[[4]][1:3000,2],1),col="darkred",type="l")


nucBottom <- c(colMeans(nuc0hrange5[[1]][rownames(nuc0hrange5[[3]]) %in% bottom1000byRNAseq,]),
            colMeans(nuc0hrange5[[2]][rownames(nuc0hrange5[[3]]) %in% bottom1000byRNAseq,]),
            colMeans(nuc0hrange5[[3]][rownames(nuc0hrange5[[3]]) %in% bottom1000byRNAseq,])
)

nuc4Bottom <- c(colMeans(nuc0hrange4[[1]][rownames(nuc0hrange4[[3]]) %in% bottom1000byRNAseq,]),
               colMeans(nuc0hrange4[[2]][rownames(nuc0hrange4[[3]]) %in% bottom1000byRNAseq,]),
               colMeans(nuc0hrange4[[3]][rownames(nuc0hrange4[[3]]) %in% bottom1000byRNAseq,])
)
plot(nuc0hrange5[[4]][1:3000,1],nucBottom[1:3000],type="l",lwd=2)
plot(nuc0hrange5[[4]][1:3000,1],runmean(nuc0hrange5[[4]][1:3000,2],1),col="darkred",type="l")
