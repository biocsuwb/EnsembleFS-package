#' Monte Carlo Feature Selection And Interdependency Discovery
#'
#' @details
#' Performs Monte Carlo Feature Selection (MCFS-ID) on a given data set.
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param ... unused arguments
#' @return A \code{\link{data.frame}} with selected features
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' fs.mcfs(data, decisions, ...)
#' }
#'
#' @import rmcfs
#'
#' @export
fs.mcfs <- function(x, y, ...){
  if (!is.data.frame(x)) x = as.data.frame(x)
  old.names <- colnames(x)
  colnames(x) <- 1:ncol(x)
  x$class <- y
  result <- mcfs(class ~., x,
                 projections = 'auto',
                 projectionSize = 'auto',
                 featureFreq = 150,
                 splits = 5,
                 splitSetSize = 1000,
                 balance = 'auto',
                 cutoffMethod = c("mean", "kmeans"),
                 cutoffPermutations = 20,
                 buildID = TRUE,
                 finalRuleset = TRUE,
                 finalCV = TRUE,
                 finalCVSetSize = 1000,
                 seed = NA,
                 threadsNumber = 5)
  result.data <- result$RI
  var.name <- names(result$data)[-length(names(result$data))]
  var.imp <- result.data[result.data$attribute %in% c(var.name),c('attribute', 'RI')]
  colnames(var.imp)[1] <- 'name'
  var.imp$name <- old.names[c(var.imp$name)]
  return(var.imp)
}
