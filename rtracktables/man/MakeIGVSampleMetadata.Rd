\name{MakeIGVSampleMetadata}
\alias{MakeIGVSampleMetadata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeIGVSampleMetadata(sampleMetadata, SampleSheet, igvdirectory)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleMetadata}{
%%     ~~Describe \code{sampleMetadata} here~~
}
  \item{SampleSheet}{
%%     ~~Describe \code{SampleSheet} here~~
}
  \item{igvdirectory}{
%%     ~~Describe \code{igvdirectory} here~~
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
function (sampleMetadata, SampleSheet, igvdirectory) 
{
    write.table("#sampleTable", file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, sep = "\t")
    colnames(sampleMetadata)[1] <- "Linking_id"
    print(colnames(sampleMetadata))
    write.table(sampleMetadata, file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = T, quote = F, append = T, 
        sep = "\t")
    BamMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[, 
        "bam"]), "SampleName"], "Bam", sep = "_"), SampleSheet[, 
        "SampleName"])
    BigWigMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[, 
        "bigwig"]), "SampleName"], "Bigwig", sep = "_"), SampleSheet[, 
        "SampleName"])
    IntervalMappings <- cbind(paste(SampleSheet[!is.na(SampleSheet[, 
        "interval"]), "SampleName"], "Interval", sep = "_"), 
        SampleSheet[, "SampleName"])
    write.table("\n#sampleMapping", file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table("#Bams", file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table(BamMappings, file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table("\n#BigWigs", file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table(BigWigMappings, file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table("\n#Intervals", file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
    write.table(IntervalMappings, file.path(igvdirectory, "SampleMetadata.txt"), 
        row.names = F, col.names = F, quote = F, append = T, 
        sep = "\t")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
