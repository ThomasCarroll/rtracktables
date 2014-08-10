#' The ChIPprofile class
#' @export
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(params="list"
         ))
