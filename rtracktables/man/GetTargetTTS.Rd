\name{GetTargetTTS}
\alias{GetTargetTTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetTargetTTS(RegionRanges, distance, distanceIn = 0, min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RegionRanges}{
%%     ~~Describe \code{RegionRanges} here~~
}
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{distanceIn}{
%%     ~~Describe \code{distanceIn} here~~
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
function (RegionRanges, distance, distanceIn = 0, min) 
{
    RegionRanges <- RegionRanges[width(RegionRanges) > min]
    TTSs <- GetAroundTTS(RegionRanges, distance, distanceIn = distanceIn)
    return(TTSs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
