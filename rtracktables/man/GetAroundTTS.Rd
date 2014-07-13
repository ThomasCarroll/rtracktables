\name{GetAroundTTS}
\alias{GetAroundTTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetAroundTTS(FullGeneBounds, distance, distanceIn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FullGeneBounds}{
%%     ~~Describe \code{FullGeneBounds} here~~
}
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{distanceIn}{
%%     ~~Describe \code{distanceIn} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (FullGeneBounds, distance, distanceIn) 
{
    require(GenomicRanges)
    PosGenes <- (FullGeneBounds[strand(FullGeneBounds) == "+"])
    NegGenes <- (FullGeneBounds[strand(FullGeneBounds) == "-"])
    TempFramePos <- cbind(as.data.frame(elementMetadata(PosGenes)), 
        width(PosGenes))
    TempFrameNeg <- cbind(as.data.frame(elementMetadata(NegGenes)), 
        width(NegGenes))
    colnames(TempFramePos)[3] <- "WidthOfOriginalGene"
    colnames(TempFrameNeg)[3] <- "WidthOfOriginalGene"
    NewPosGenes <- GRanges(seqnames = seqnames(PosGenes), IRanges(start = (as.vector(end(ranges(PosGenes)))) - 
        distanceIn, end = (as.vector(end(ranges(PosGenes)))) + 
        distance), strand = strand(PosGenes), TempFramePos)
    NewNegGenes <- GRanges(seqnames = seqnames(NegGenes), IRanges(start = (as.vector(start(ranges(NegGenes))) - 
        distance), end = (as.vector(start(ranges(NegGenes)))) + 
        distanceIn), strand = strand(NegGenes), TempFrameNeg)
    names(NewPosGenes) <- names(PosGenes)
    names(NewNegGenes) <- names(NegGenes)
    AllPromoters <- c(NewPosGenes, NewNegGenes)
    return(AllPromoters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
