#' Run end-to-end Ensemble for comparison of feature selection methods.
#'
#' @details
#' Ensemble for comparison of feature selection methods dedicated to high-throughput
#' sequencing data.
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param methods A \code{\link{vector}}  with feature selection methods available in this library for comparison
#' @param method.cv validation method \code{kfoldcv} for cross-validation \code{k-fold} or \code{rsampling} for \code{random sampling}
#' @param params.cv A \code{\link{list}} with the following fields:
#' \itemize{
#'   \item \code{k} --  the number of groups that a given data sample is to be split into, not less than 3
#'   \item \code{test.size} -- testing set size for random sampling validation
#'   \item \code{iter} -- the number of validation repetitions
#' }
#' @param level.cor cutoff level of correlated variables. If equal to 1 is not performed
#' @param params A \code{\link{list}} with the following fields:
#' \itemize{
#'   \item \code{adjust} -- method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"} is recommended for FDR, see Details) for MDFS1D, MDFS2D and U-test
#'   \item \code{feature.number} -- number of attributes to select. Must not exceed \code{ncol(x)}
#'   \item \code{alpha} -- significance level for MDFS1D, MDFS2D and U-test
#'   \item \code{use.cuda} -- whether to use CUDA acceleration (must be compiled with CUDA) for MDFS2D method
#'   \item \code{cutoff.method} -- cutoff method MCFS: "permutations", "criticalAngle", "kmeans", "mean", "contrast"
#' }
#' @param asm A \code{\link{vector}} with enumeration method for which to calculate Lustgarten’s stability measure
#' @param model A \code{\link{vector}} with enumeration method for which to training and testing model \code{RandomForest}
#' @return
#' \itemize{
#'   \item \code{selected.feature} -- A \code{\link{list}} with the result of feature selection for the selected feature selection method
#'   \item \code{ranking.feature} -- A \code{\link{list}} with the result of the rating of the variables that were most often performed in each iteration of cross-validation
#'   \item \code{stability} -- A \code{\link{data.frame}} with the result of stability of selection of feature for the selected selection method
#'   \item \code{model} -- A \code{\link{data.frame}} with the result of constructing a random forest model for the selected feature selection method
#' }
#'
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' ensembleFS(data,
#'            decisions,
#'            methods = c('fs.utest', 'fs.mrmr'),
#'            method.cv = 'kfoldcv',
#'            params.cv = list(k = 3, iter = 10),
#'            level.cor = 0.75,
#'            params = list(adjust = 'fdr', feature.number = 10, alpha = 0.05),
#'            asm = c('fs.utest', 'fs.mrmr'),
#'            model = c('fs.utest', 'fs.mrmr')
#'            )
#'
#' }
#'
#' @importFrom stats var
#' @export
ensembleFS <- function(x,
                        y,
                        methods = c('fs.utest'),
                        method.cv = 'kfoldcv',
                        params.cv = list(k = 3, niter = 5),
                        level.cor = 1,
                        params = list(adjust = 'holm',
                                      feature.number = 10,
                                      alpha = 0.05,
                                      use.cuda = FALSE,
                                      cutoff.method = c("kmeans")),
                        asm = c('fs.utest'),
                        model = c('fs.utest')){

  if(!is.data.frame(x)){x <- as.data.frame(x)}

  if (length(y) != nrow(x)) {
    stop('Length of decision is not equal to the number of rows in data.')
  }

  if (!any(y == 0)) {
    y <- y - as.integer(1)
  }

  if (!all(y == 0 | y == 1)) {
    stop('Decision must be binary.')
  }

  if (all(y == 0) || all(y == 1)) {
    stop('Both classes have to be represented.')
  }

  nums <- unlist(lapply(x, is.numeric))
  if(FALSE %in% nums){
    stop('Columns in data must be of type numeric!')
  }

  for(i in asm){
    if(!(is.element(i, methods))){
      stop('methods from the asm argument must be listed in methods argument')
    }
  }

  for(i in model){
    if(!(is.element(i, methods))){
      stop('methods from the model argument must be listed in methods argument')
    }
  }

  x = x[, !apply(x, 2, var) == 0]

  list.index.cross <- cross.val(x, y, method.cv, params.cv)
  feature.selection.result <- list()
  for(method in methods){
    result <- feature.selection.cv(x, y, method, list.index.cross, params = params)
    feature.selection.result <- append(feature.selection.result, list(result))
  }
  names(feature.selection.result) <- methods
  if(level.cor != 1){
    result.uncor.var <- list()
    for(i in 1:length(feature.selection.result)){
      result <- corelation.removed(x, feature.selection.result[[i]], list.index.cross, level.cor)
      result.uncor.var <- append(result.uncor.var, list(result))
    }
    names(result.uncor.var) <- methods
    feature.selection.result <- result.uncor.var
  }

  result.ranking <- data.frame()
  for(i in 1:length(feature.selection.result)){
    result <- ranking.feature(feature.selection.result[[i]])
    result$method <- substring(names(feature.selection.result[i]), 4)
    result.ranking <- rbind(result.ranking, result)
  }

  stability.selection.result <- data.frame()
  for(i in 1:length(feature.selection.result)){
    if(names(feature.selection.result[i]) %in% asm){
      result <- stability.selection.top.var(feature.selection.result[[i]], list.index.cross)
      result$method <- substring(names(feature.selection.result[i]), 4)
      stability.selection.result <- rbind(stability.selection.result, result)
    }
  }

  metrics.model <- data.frame()
  for(i in 1:length(feature.selection.result)){
    if(names(feature.selection.result[i]) %in% model){
      result <- model.result.top.var(x, y,  feature.selection.result[[i]], list.index.cross)
      result$method <- substring(names(feature.selection.result[i]), 4)
      metrics.model <- rbind(metrics.model, result)
    }
  }
  return(out = list(selected.feature = feature.selection.result,
                    ranking.feature = result.ranking,
                    stability = stability.selection.result,
                    model = metrics.model))
}
