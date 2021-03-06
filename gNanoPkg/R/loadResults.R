#' Load results file
#'
#' Provides a simple way to load results files and turns
#' them into a data.frame that can be manipulated
#'
#' @param model
#' @param resultsRoot
#' @param summary
#'
#' @return
#' @export
#'
#' @examples
loadResults = function(model = c(paste0("g-",1:7),
                                 paste0("ln-",1:7),
                                 paste0("sn-",1:7)),
                       resultsRoot = "../systematic/",
                       summary = FALSE){

  model = match.arg(model)

  resultsFile = glue("{resultsRoot}{model}.Rda")
  results = load(resultsFile)

  if(summary){
    return(simSummary)
  }

  return(as.data.frame(sim.sample[[1]]))
}
