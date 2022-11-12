#' Minimum Redundancy Maximal Relevancy
#'
#' @details
#' Performs MRMR for feature selection
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param params number of attributes to select. Must not exceed, feature.number there cannot be more than the number of columns
#' @return A \code{\link{data.frame}} with selected features
#'
#' @examples
#'
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' fs.mrmr(data, decisions, params = list(feature.number = 100))
#' }
#'
#' @import mRMRe
#'
#' @export
fs.mrmr <- function(x, y, params = list(feature.number = 100)){
  if (!is.data.frame(x)) x = as.data.frame(x)
  if (params$feature.number > ncol(x)){
    stop('fs.mrmr :feature.number there cannot be more than the number of columns')
  }
  x <- sapply(x, as.double)
  x <- cbind(class = y, x)
  old.names.x = colnames(x)
  x = transform(x, class = as.numeric(class))
  new.names.x = as.character(sprintf("a%d", 1:ncol(x)))
  names(x) = new.names.x
  x <- as.data.frame(x)
  feature.data <- mRMR.data(data = x)
  results <- mRMR.ensemble(data = feature.data,
                           target_indices = 1,
                           solution_count = 1,
                           feature_count = params$feature.number)
  var.imp.new.names = apply(solutions(results)[[1]], 2,
                            function(x, y) { return(y[x]) },
                            y = featureNames(feature.data))
  scores = as.data.frame(unlist(scores(results)))
  name = old.names.x[match(var.imp.new.names,names(x))]
  var.imp = cbind(name, scores)
  names(var.imp) <- c('name', 'score')
  return(var.imp)
}
