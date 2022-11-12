#' Build MultiDimensional Feature Selector from IGs uses GPU,  is a parallel computing platform CUDA
#'
#' @details
#' Build MultiDimensional Feature Selector from IGs, requires mandatory installation of CUDA
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param params A \code{\link{list}} with the following fields:
#' \itemize{
#'   \item \code{adjust} -- method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details)
#'   \item \code{alpha} -- significance level for MDFS1D, MDFS2D and U-test
#'   \item \code{use.cuda} -- whether to use CUDA acceleration (must be compiled with CUDA) for MDFS2D method
#' }
#' @return A \code{\link{data.frame}} with selected features and p.value
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' fs.mdfs.2D(data, decisions, params = list(adjust = 'holm', use.cuda = FALSE))
#' }
#'
#' @import MDFS
#' @importFrom stats p.adjust
#'
#' @export
fs.mdfs.2D <- function(x, y, params = list(adjust = 'holm', alpha = 0.05, use.cuda = FALSE)){
  if (!is.data.frame(x)) data = as.data.frame(x)
  dim0 = 2
  div0 = 1
  adjust = params$adjust
  use.cuda = params$use.cuda
  alpha = params$alpha
  result = MDFS(data = x, decision = y, dimensions = dim0, divisions = div0, use.CUDA = use.cuda,
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

