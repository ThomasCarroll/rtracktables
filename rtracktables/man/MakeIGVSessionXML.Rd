\name{MakeIGVSessionXML}
\alias{MakeIGVSessionXML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeIGVSessionXML(files, fileType, names, igvdirectory, XMLname, genomeName, locusName = "All")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
%%     ~~Describe \code{files} here~~
}
  \item{fileType}{
%%     ~~Describe \code{fileType} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{igvdirectory}{
%%     ~~Describe \code{igvdirectory} here~~
}
  \item{XMLname}{
%%     ~~Describe \code{XMLname} here~~
}
  \item{genomeName}{
%%     ~~Describe \code{genomeName} here~~
}
  \item{locusName}{
%%     ~~Describe \code{locusName} here~~
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
function (files, fileType, names, igvdirectory, XMLname, genomeName, 
    locusName = "All") 
{
    i <- 1
    require(XML)
    files <- gsub(".xls", ".bed", files)
    Output <- file.path(igvdirectory, paste(XMLname, ".xml", 
        sep = ""))
    resources <- vector("list", length = length(files))
    if (fileType == "Bam") {
        NewName <- paste(names[i], "_Bam", sep = "")
    }
    if (fileType == "Bigwig") {
        NewName <- paste(names[i], "_Bigwig", sep = "")
    }
    if (fileType == "Macs") {
        NewName <- paste(names[i], "_Macs", sep = "")
    }
    if (fileType == "MacsDiff") {
        NewName <- paste(names[i], "_MacsDiff", sep = "")
    }
    if (fileType == "TPICs") {
        NewName <- paste(names[i], "_TPICs", sep = "")
    }
    if (fileType == "Sicer") {
        NewName <- paste(names[i], "_Sicer", sep = "")
    }
    GlobalNode <- newXMLNode("Global", attrs = c(genome.value = genomeName, 
        groupTracksBy = "Linking_id", locus = locusName, version = 3))
    ResourcesNode <- newXMLNode("Resources", parent = GlobalNode)
    for (i in 1:length(resources)) {
        resources[[i]] <- newXMLNode("Resource", parent = ResourcesNode, 
            attrs = c(label = NewName[i], name = NewName[i], 
                path = relativePath(files[i], Output), relativePath = T))
    }
    MetaDataNode <- newXMLNode("Resource", parent = ResourcesNode, 
        attrs = c(name = "SampleMetadata", path = relativePath(file.path(igvdirectory, 
            "SampleMetadata.txt"), Output), relativePath = T))
    PanelDataNode <- newXMLNode("Panel", attrs = c(height = "350", 
        name = "DataPanel", width = "1115"), parent = GlobalNode)
    if (fileType == "Bam") {
        TrackNode <- newXMLNode("Track", attrs = c(altColor = "0,0,178", 
            color = "0,0,178", colorOption = "UNEXPECTED_PAIR", 
            displayMode = "EXPANDED", featureVisibilityWindow = "-1", 
            fontSize = "10", id = relativePath(files[i], Output), 
            name = NewName[i], showDataRange = "true", sortByTag = "", 
            visible = "true"), parent = PanelDataNode)
    }
    if (fileType == "Bigwig") {
        TrackNode <- newXMLNode("Track", attrs = c(altColor = "0,0,178", 
            autoscale = "true", color = "0,0,178", displayMode = "COLLAPSED", 
            featureVisibilityWindow = "-1", fontSize = "10", 
            id = relativePath(files[i], Output), name = NewName[i], 
            renderer = "BAR_CHART", showDataRange = "true", visible = "true", 
            windowFunction = "mean"), parent = PanelDataNode)
        DisplayRangeNode <- newXMLNode("DataRange", attrs = c(baseline = "0.0", 
            drawBaseline = "true", flipAxis = "false", maximum = "50", 
            minimum = "5", type = "LINEAR"), parent = TrackNode)
    }
    if (fileType == "Macs" | fileType == "TPICs" | fileType == 
        "Sicer" | fileType == "MacsDiff") {
        TrackNode <- newXMLNode("Track", attrs = c(altColor = "0,0,178", 
            color = "0,0,178", displayMode = "COLLAPSED", featureVisibilityWindow = "-1", 
            fontSize = "10", height = "45", id = relativePath(files[i], 
                Output), name = NewName[i], renderer = "BASIC_FEATURE", 
            showDataRange = "true", sortable = "false", visible = "true", 
            windowFunction = "count"), parent = PanelDataNode)
    }
    saveXML(GlobalNode, file = Output)
    return(Output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
