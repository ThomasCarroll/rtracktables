#' The ChIPprofile class
#'
#' @slot profile list
setClass("ChIPprofile",contains = "SummarizedExperiment",
         slots=c(profile="list"
         ))
