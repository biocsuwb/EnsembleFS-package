#' List methods
#'
#' @details
#' List methods available in the library
#'
#' @export
list.methods <- function(){
  everything <- sort(getNamespaceExports("ensembleFS"))
  message("All feature selection algorithm wrappers in benchmarkFS:\n")
  print(everything[grepl(pattern="^[f]s", everything)])
}

