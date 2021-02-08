#' Showing result ensemble
#'
#' @details
#' The function creates a graph showing Lustgarten's stability measure ASM (N)
#' for top-N variables N = 5,10,15,20,30,40,50,75,100 for selected FS
#' The function creates a graph showing the dependence of the selected metrics model(N)
#' N = 5,10,15,20,30,40,50,75,100 for selected FS
#'
#'
#' @param x  A \code{\link{data.frame}} with data result metrics model and  Lustgarten’s stability measure
#' @param y what is showing , please selected : stability, acc, auc, mcc.
#' @return graph showing
#'
#' @examples
#'
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' result <- ensembleFS(data,
#'            decisions,
#'            methods = c('fs.utest', 'fs.mrmr'),
#'            method = 'kfoldcv',
#'            params.cv = list(k = 3, iter = 10),
#'            level.cor = 0.75,
#'            params = list(adjust = 'SGoF', feature.number = 10, alpha = 0.05),
#'            asm = c('fs.utest', 'fs.mrmr'),
#'            model = c('fs.utest', 'fs.mrmr')
#'            )
#'
#'
#'
#' graph.result(result$stability, 'stability')
#'
#' graph.result(result$model, 'acc')
#'
#' }
#' @rdname graph.result
#'
#' @import ggplot2
#' @export
graph.result <- function(x, y){
  if(!(y %in% c('stability', 'acc', 'auc', 'mcc'))){
    stop('available stability, acc, auc, mcc')
  }
  if(y == 'stability'){
    result <- ggplot(data = x, aes(x = N, y = stability.asm, group= method ,color = method)) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks= c(seq(0,100, by = 10))) +
      labs(title= "ASM similarity measure between feature subsets vs top N variables.", y="ASM", x = "N")
    return(result)
  }
  if(y == 'acc'){
    result <- ggplot(data = x, aes(x = N, y = mean.acc, group= method ,color = method)) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks= c(seq(0,100, by = 10))) +
      labs(title= "The accuracy vs top N variables.", y="ACC", x = "N")

    return(result)
  }
  if(y == 'auc'){
    result <- ggplot(data = x, aes(x = N, y = mean.auc, group= method ,color = method)) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks= c(seq(0,100, by = 10))) +
      labs(title= "Area under the ROC curve vs top N variables.", y="AUC", x = "N")

    return(result)
  }
  if(y == 'mcc'){
    result <- ggplot(data = x, aes(x = N, y = mean.mcc, group= method ,color = method)) +
      geom_line() +
      geom_point() +
      scale_x_continuous(breaks= c(seq(0,100, by = 10))) +
      labs(title= "Matthews Correlation Coefficient vs top N variables.", y="MCC", x = "N")

    return(result)
  }
}
