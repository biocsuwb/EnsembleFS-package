#' Build MultiDimensional Feature Selector from IGs
#'
#' @details
#' Build MultiDimensional Feature Selector from IGs
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param params method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"}
#' is recommended for FDR, see Details)
#' @return A \code{\link{data.frame}} with selected features and p.value
#'
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' fs.mdfs.1D(data, decisions, params = list(adjust = 'holm', alpha = 0.05))
#' }
#'
#'
#' @import MDFS
#' @importFrom stats p.adjust
#' @export
fs.mdfs.1D <- function(x, y, params = list(adjust = 'holm', alpha = 0.05)){
  if (!is.data.frame(x)) data = as.data.frame(x)
  dim0 = 1
  div0 = 3
  adjust = params$adjust
  alpha = params$alpha 
  result = MDFS(data = x, decision = y, dimensions = dim0, divisions = div0, use.CUDA = FALSE,
                p.adjust.method = adjust)
  var.names = names(x)
  index.imp = RelevantVariables(result$MDFS,
                                level = alpha,
                                p.adjust.method = adjust)
  var.imp.frame = data.frame(name = var.names,
                             Pvalue = result$p.value,
                             adjustPvalue = result$adjusted.p.value)[index.imp,]
  var.imp = var.imp.frame[order(var.imp.frame$adjustPvalue, decreasing = F),]
  return(var.imp)
}

