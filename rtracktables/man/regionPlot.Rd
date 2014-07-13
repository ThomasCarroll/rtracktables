\name{regionPlot}
\alias{regionPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regionPlot(bamFile, testRanges, nOfWindows = 100, FragmentLength = 150, style = "point", distanceAround = 1500, distanceInRegionStart = 1500, distanceOutRegionStart = 1500, distanceInRegionEnd = 1500, distanceOutRegionEnd = 1500, paired = F, normalize = "RPM", plotBy = "coverage", removeDup = F, verbose = T, format = "bam", seqlengths = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFile}{
%%     ~~Describe \code{bamFile} here~~
}
  \item{testRanges}{
%%     ~~Describe \code{testRanges} here~~
}
  \item{nOfWindows}{
%%     ~~Describe \code{nOfWindows} here~~
}
  \item{FragmentLength}{
%%     ~~Describe \code{FragmentLength} here~~
}
  \item{style}{
%%     ~~Describe \code{style} here~~
}
  \item{distanceAround}{
%%     ~~Describe \code{distanceAround} here~~
}
  \item{distanceInRegionStart}{
%%     ~~Describe \code{distanceInRegionStart} here~~
}
  \item{distanceOutRegionStart}{
%%     ~~Describe \code{distanceOutRegionStart} here~~
}
  \item{distanceInRegionEnd}{
%%     ~~Describe \code{distanceInRegionEnd} here~~
}
  \item{distanceOutRegionEnd}{
%%     ~~Describe \code{distanceOutRegionEnd} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{plotBy}{
%%     ~~Describe \code{plotBy} here~~
}
  \item{removeDup}{
%%     ~~Describe \code{removeDup} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{format}{
%%     ~~Describe \code{format} here~~
}
  \item{seqlengths}{
%%     ~~Describe \code{seqlengths} here~~
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
function (bamFile, testRanges, nOfWindows = 100, FragmentLength = 150, 
    style = "point", distanceAround = 1500, distanceInRegionStart = 1500, 
    distanceOutRegionStart = 1500, distanceInRegionEnd = 1500, 
    distanceOutRegionEnd = 1500, paired = F, normalize = "RPM", 
    plotBy = "coverage", removeDup = F, verbose = T, format = "bam", 
    seqlengths = NULL) 
{
    if (!verbose) {
        suppressMessages(runRegionPlot())
    }
    result <- runRegionPlot(bamFile, testRanges, nOfWindows, 
        FragmentLength, style, distanceAround, distanceInRegionStart, 
        distanceOutRegionStart, distanceInRegionEnd, distanceOutRegionEnd, 
        paired, normalize, plotBy, removeDup, format, seqlengths)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
