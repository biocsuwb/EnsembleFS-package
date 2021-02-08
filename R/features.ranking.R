#' Function for ranking selected variables
#'
#' @details
#' Ranking is calculated based on selected variables that are more common in cross-validation iterations.
#'
#' @param list.selected.var A \code{\link{list}} with selected variables in cross-validation
#' @return A \code{\link{data.frame}} with the result of the rating of the variables that were most often performed in each iteration of cross-validation
#'
#' @examples
#'
#' \dontrun{
#'
#' class <- data$class
#' data$class <- NULL
#'
#' list.index.cross <- cross.val(x = data,
#'                              y = decisions,
#'                              method = 'kfoldcv',
#'                              params.cv = list(niter = 10, k = 3))
#'
#'
#' list.selected.var <- feature.selection(x = data,
#'                              y = class,
#'                              method = 'fs.utest',
#'                              list.index.cross = list.index.cross,
#'                              params = list(adjust = 'holm'))
#'
#' ranking.var <- ranking.feature(list.selected.var = list.selected.var)
#'
#' }
#'
#' @export
ranking.feature <- function(list.selected.var){
  var.list <- c()
  for(i in 1:length(list.selected.var)){
    var.list <- append(var.list, as.character(list.selected.var[[i]]$name))
  }
  var.info <- as.data.frame(table(var.list))
  colnames(var.info) <- c('biomarker.name', 'frequency')
  result.ranking.var <- var.info[order(var.info$frequency, decreasing=TRUE),]
  return(result.ranking.var)
}
