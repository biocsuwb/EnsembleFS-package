#' Test U Manna-Whitneya (U-test) for feature selection
#'
#' @details
#' In case of FDR control it is recommended to use Benjamini-Hochberg-Yekutieli
#' p-value adjustment
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param params method as accepted by \code{\link[stats]{p.adjust}} (\code{"BY"}
#' is recommended for FDR, see Details)
#' @return A \code{\link{data.frame}} with selected features and p.value
#'
#' @examples
#'
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' fs.utest(data, decisions, params = list(adjust = 'holm', alpha = 0.05))
#' }
#'
#' @importFrom stats wilcox.test p.adjust
#' 
#' @export
fs.utest <- function(x, y, params = list(adjust = 'holm', alpha = 0.05)){
  if(!is.data.frame(x)){ x <- as.data.frame(x)}
  x$class <- y
  data.class1 <- x[x$class == 0,1:ncol(x)]
  data.class2 <- x[x$class == 1,1:ncol(x)]
  data.class1$class <- NULL
  data.class2$class <- NULL
  doUtest = function(i,x,y){
    data.part = wilcox.test(data.class1[,i],
                            data.class2[,i],
                            alternative = "two.sided",
                            paired = FALSE)
    p.data.part=data.part$p.value
    return(p.data.part)
  }
  ncol.class1 = ncol(data.class1)
  p.value = sapply(1:ncol.class1, doUtest, x=data.class1, y=data.class2)
  var.p.value = as.data.frame(cbind(colnames(data.class1),p.value))
  names(var.p.value) = c("name","Pvalue")
  var.p.value$Pvalue = as.numeric(as.character((var.p.value$Pvalue)))
  var.p.value.sort = var.p.value[order(var.p.value$Pvalue, decreasing=F),]
  adjust = params$adjust
  adjust.p.value = p.adjust(var.p.value.sort$Pvalue, method = adjust)
  var.imp.all = as.data.frame(cbind(var.p.value.sort, adjust.p.value))
  names(var.imp.all) = c("name","Pvalue","adjustPvalue")
  row.names(var.imp.all) = NULL
  var.imp = var.imp.all[which(var.imp.all$adjustPvalue < params$alpha),]
  return(var.imp)
}
