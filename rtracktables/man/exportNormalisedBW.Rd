\name{exportNormalisedBW}
\alias{exportNormalisedBW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exportNormalisedBW(bamFile, qc, normaliseTo = "blacklisted")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFile}{
%%     ~~Describe \code{bamFile} here~~
}
  \item{qc}{
%%     ~~Describe \code{qc} here~~
}
  \item{normaliseTo}{
%%     ~~Describe \code{normaliseTo} here~~
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
function (bamFile, qc, normaliseTo = "blacklisted") 
{
    require(GenomicAlignments)
    require(rtracklayer)
    extendBy <- fragmentlength(qc)
    message("Reading tags from ", bamFile, appendLF = FALSE)
    if (normaliseTo == "blacklisted") {
        totalReads <- qc@FlagAndTagCounts["MapQPass"] - qc@CountsInFeatures$BlackList
    }
    if (normaliseTo == "Total") {
        totalReads <- qc@FlagAndTagCounts["Mapped"]
    }
    if (normaliseTo == "UniqueTotal") {
        totalReads <- qc@FlagAndTagCounts["Mapped"] - qc@FlagAndTagCounts["Duplicates"]
    }
    total <- readGAlignmentsFromBam(bamFile)
    message("..done")
    message("Read in ", length(total), " reads")
    message("Extending reads to fragmentlength of ", extendBy, 
        " ..", appendLF = FALSE)
    temp <- resize(as(total, "GRanges"), extendBy, "start")
    message("..done")
    rm(total)
    gc()
    message("Calculating coverage..", appendLF = FALSE)
    genomeCov <- coverage(temp)
    rm(temp)
    message("..done")
    message("Normalised coverage..", appendLF = FALSE)
    genomeCov <- (genomeCov/totalReads) * 1e+06
    message("..done")
    message("Exporting coverage..", appendLF = FALSE)
    export.bw(genomeCov, file.path(dirname(dirname(bamFile)), 
        gsub("\\.bam", "Normalised\\.bw", basename(bamFile))))
    message("..done")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
