\name{GetGene}
\alias{GetGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetGene(FullGeneBounds, distanceFromStart, distanceFromEnd, min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FullGeneBounds}{
%%     ~~Describe \code{FullGeneBounds} here~~
}
  \item{distanceFromStart}{
%%     ~~Describe \code{distanceFromStart} here~~
}
  \item{distanceFromEnd}{
%%     ~~Describe \code{distanceFromEnd} here~~
}
  \item{min}{
%%     ~~Describe \code{min} here~~
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
function (FullGeneBounds, distanceFromStart, distanceFromEnd, 
    min) 
{
    require(GenomicRanges)
    FullGeneBounds <- FullGeneBounds[width(FullGeneBounds) > 
        min]
    PosGenes <- (FullGeneBounds[strand(FullGeneBounds) == "+"])
    NegGenes <- (FullGeneBounds[strand(FullGeneBounds) == "-"])
    NewPosGenes <- GRanges(seqnames = seqnames(PosGenes), IRanges(start = (as.vector(start(ranges(PosGenes)))) - 
        distanceFromStart, end = (as.vector(end(ranges(PosGenes)))) + 
        distanceFromEnd), strand = strand(PosGenes), elementMetadata(PosGenes))
    NewNegGenes <- GRanges(seqnames = seqnames(NegGenes), IRanges(start = (as.vector(start(ranges(NegGenes))) - 
        distanceFromEnd), end = (as.vector(end(ranges(NegGenes)))) + 
        distanceFromStart), strand = strand(NegGenes), elementMetadata(NegGenes))
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
