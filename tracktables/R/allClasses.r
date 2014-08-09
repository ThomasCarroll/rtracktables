#' The ChIPprofile class
#' @export
#' @slot profile list
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(profile="list"
         ))
