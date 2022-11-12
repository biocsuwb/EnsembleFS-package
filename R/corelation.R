#' @importFrom stats cor
corelation <- function(data, level.cor, start.var){
  .cor = function(m, data, level.cor, start.var){
    if (!is.data.frame(data)) data = as.data.frame(data)
    core = abs(cor(as.numeric(data[,start.var]), as.numeric(data[,m]), method = "spearman"))
    if(core > level.cor  & m != start.var){
      return(m)
    }
  }
  cor.list = lapply(start.var:ncol(data), function(m) .cor(m, data, level.cor, start.var))
  return(cor.list)
}


#' Function to search highly correlated variables
#'
#' @details
#' Searching highly  correlated variables
#'
#' @param data input data where columns are variables and rows are observations (all numeric)
#' @param level.cor cutoff level of correlated variables. If equal to 1 is not performed
#' @return A \code{\link{vector}} with indexes of highly correlated variables
#'
#' @export
corelation.search <- function(data, level.cor){
  listcor <- lapply(1:ncol(data), function(start) corelation(data, level.cor, start.var = start))
  return(unique(unlist(listcor)))
}


#' Function to remove highly correlated variables
#'
#' @details
#' Removes highly correlated variables after selection performed in cross-validation
#'
#' @param data input data where columns are variables and rows are observations (all numeric)
#' @param list.selected.var A \code{\link{list}} with selected variables in cross-validation
#' @param list.index.cross A \code{\link{list}} with indexes obtained in cross-validation
#' @param level.cor cutoff level of correlated variables. If equal to 1 is not performed
#' @return A \code{\link{list}} with selected features for each iteration of cross-validation
#'
#' @examples
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
#' list.selected.var <- feature.selection(x = data,
#'                                        y = class,
#'                                        method = 'fs.utest',
#'                                        list.index.cross = indexes,
#'                                        params = list(adjust = 'fdr'))
#'                              
#' list.selected.var.uncor <- corelation.removed(data,
#'                                               list.selected.var,
#'                                               list.index.cross,
#'                                               0.75)
#'
#' }
#'
#' @export
corelation.removed <- function(data, list.selected.var, list.index.cross, level.cor){
  if (!is.data.frame(data)) data = as.data.frame(data)
  niter = length(list.index.cross)
  ncross = length(list.index.cross[[1]]$training)
  cor.rm <- function(m, data, list.selected.var, list.index.cross, niter, ncross, level.cor){
    index = expand.grid(j=1:ncross, p=1:niter)
    p <- index[m,2]
    j <- index[m,1]
    idx.cross <- list.index.cross[[p]]
    data.cross <- na.omit(data[idx.cross$training[[j]],])
    var <- list.selected.var[[m]]
    if(nrow(var) == 0 | nrow(var) == 1){
      return(var)
    }
    else{
      idx.cor.var <- corelation.search(data.cross[,c(var$name)], level.cor)
     }
    if(is.null(idx.cor.var)){
      return(var)
    }
    else{
      var <- var[-(idx.cor.var),]
      return(var)
    }
  }
  N = niter*ncross
  var.uncor.iter = lapply(1:N, function(m) cor.rm(m,
                                                  data,
                                                  list.selected.var,
                                                  list.index.cross,
                                                  niter,
                                                  ncross,
                                                  level.cor))
}

