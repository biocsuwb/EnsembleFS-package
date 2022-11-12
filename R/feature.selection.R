#' Selection of features in cross-validation
#'
#' @details Uses the methods of selection of features  available in this library
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param method the name of the selection method for traits available in this library
#' @param list.index.cross A \code{\link{list}} with indexes obtained in cross-validation
#' @param params A \code{\link{list}} parameters required for this method of feature selection
#' @return  A \code{\link{list}} with selected features for each iteration of cross-validation
#' @examples
#'
#' \dontrun{
#'
#' class <- data$class
#' data$class <- NULL
#'
#' list.index.cross <- cross.val(x = data,
#'                               y = class,
#'                               method= 'kfoldcv',
#'                               params.cv = list(k = 3, niter = 10))
#'
#' list.selected.var <- feature.selection.cv(x = data,
#'                              y = class,
#'                              method = 'fs.utest',
#'                              list.index.cross = list.index.cross,
#'                              params = list(adjust = 'holm', alpha = 0.05))
#'
#' }
#' @export
feature.selection.cv = function(x, y, method, list.index.cross, params = list(adjust = 'holm', alpha = 0.05)){
  if (!is.data.frame(data)) x = as.data.frame(x)
  niter = length(list.index.cross)
  ncross = length(list.index.cross[[1]]$training)
  var.imp.selection = function(m, x, y, method, list.index.cross ,niter, ncross, params = list()){
    index = expand.grid(j=1:ncross, p=1:niter)
    p = index[m,2]
    j = index[m,1]
    list.index.cross.train = list.index.cross[[p]]
    x$class <- y
    data.train.cross = na.omit(x[list.index.cross.train$training[[j]],])
    decision.train.cross <- data.train.cross$class
    data.train.cross$class <- NULL
    params <- list(data.train.cross, decision.train.cross, params)
    var.imp <- do.call(method, params)
    return(var.imp)

  }

  N = niter*ncross
  var.each.iter = lapply(1:N, function(m) var.imp.selection(m,
                                                            x,
                                                            y,
                                                            method = method,
                                                            list.index.cross,
                                                            niter,
                                                            ncross,
                                                            params = params))
  return(var.each.iter)
}


#' Selection of features
#'
#' @details Uses the methods of selection of features  available in this library
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param methods the name of the selection method for traits available in this library
#' @param connection.method method of combining selection results if more than one method is selected, options \code{'sum'} for the sum of all selected filtration methods or \code{'intersect'} for part of the common
#' @param params A \code{\link{list}} parameters required for this method of feature selection
#' @return  A  new \code{\link{data.frame}} consisting of the selected functions and the decision variable
#' @examples
#'
#' \dontrun{
#'
#' class <- data$class
#' data$class <- NULL
#'
#' feature.selection(x = data,
#'                   y = class,
#'                   methods = c('fs.utest','fs.mdfs.1D', 'fs.mrmr', 'fs.mcfs'),
#'                   connection.method = 'sum',
#'                   params = list(feature.number=100, adjust = 'holm', alpha = 0.05))
#' }
#' @export
feature.selection <- function(x,
                              y,
                              methods,
                              connection.method = 'sum',
                              params = list(adjust='holm',
                                            feature.number=10,
                                            alpha=0.05)){
  if(!is.data.frame(x)){x <- as.data.frame(x)}

  if (!(connection.method %in% c('sum', 'intersect'))) {
    stop('Unknown  connection.method')
  }

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

  result <- list()
  params <- list(x, y, params)
  for(method in methods){
    var <- do.call(method, params)
    result <- append(result, list(as.vector(var[,'name'])))
  }
  if(connection.method == 'intersect'){
    var.imp <- Reduce(intersect, result)
  }
  else if(connection.method == 'sum'){
    var.imp <- unique(unlist(result))
  }
  result.data <- x[,var.imp]
  result.data$class <- y
  return(result.data)
}



