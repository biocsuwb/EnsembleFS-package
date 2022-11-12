#' Getting information about genes from gprogiler2
#'
#' @details
#' Getting information about genes from gprogiler2, requires internet connection
#'
#' @param rel.var  A \code{\link{vector}} with gene names
#' @return A \code{\link{data.frame}} with gene information
#'
#' @examples
#' \dontrun{
#'
#' gene <- c('TMC8')
#' get.info.gprofiler(gene)
#' }
#'
#' @import gprofiler2
#' @export
get.info.gprofiler = function(rel.var){
  df = list()
  for (i in 1:length(rel.var)){
    if(try(is.null(gost(query = rel.var[i])))){next}
    all.info = gost(query = rel.var[i])
    df.res = all.info$result
    df[[i]] = data.frame(term = rel.var[i], source = df.res$source, term.ID = df.res$term_id, term.name = df.res$term_name)
  }
  df = do.call(rbind,df)
  return(df)
}


#' Get the best genes
#'
#' @details
#' Get the best genes out obtained in \code{ensembleFS}
#'
#' @param list.imp.var.cv the main genes obtained in \code{ensembleFS}
#' @param level.freq how many times a gene has occurred in of the feature subsets
#' @param number.gene number of genes
#' @return A \code{\link{list}} with top genes for each method specified in \code{ensembleFS}
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' result <- ensembleFS(data,
#'             decisions,
#'             methods = c('fs.utest', 'fs.mrmr'),
#'             method.cv = 'kfoldcv',
#'             params.cv = list(k = 3, iter = 10),
#'             level.cor = 0.75,
#'             params = list(adjust = 'fdr', feature.number = 10, alpha = 0.05),
#'             asm = c('fs.utest', 'fs.mrmr'),
#'             model = c('fs.utest', 'fs.mrmr')
#'             )
#'
#' gene.top <- get.top.gene(result$selected.feature, 15 , 20)
#'
#' }
#'
#' @export
get.top.gene <- function(list.imp.var.cv, level.freq, number.gene){
  result <- c()
  name.method <- c()
  for(name in names(list.imp.var.cv)){
    var.one.method <- list.imp.var.cv[[name]]
    var.imp <- list()
    for(i in var.one.method){
      if(length(i$name) < number.gene) {var.imp <- append(var.imp, i$name)}
      else{var.imp <- append(var.imp, i$name[1:number.gene])}
    }
    var.frequency <- as.data.frame(table(unlist(var.imp)))
    colnames(var.frequency) <- c('name', 'frequency')
    var.frequency <- var.frequency[order(var.frequency$frequency, decreasing = TRUE),]
    var.frequency <- var.frequency[var.frequency$frequency >= level.freq,]
    name.method <- append(name.method, substring(name, 4))
    result <- append(result, list(var.frequency$name))
  }
  names(result) <- name.method
  return(result)
}


#' Getting information about top genes from gprogiler2
#'
#' @details
#' Getting information about top genes from gprofiler2 obtained from ensembleFS.
#' Requires internet connection.
#'
#' @param gene.top A \code{\link{list}} with top genes for each method specified in \code{ensembleFS}
#' @param condition.methods A way to combine the best genes for further analysis obtained by different methods in \code{ensembleFS}. For part of the common \code{intersect}, for union from methods \code{union}.
#' @return A \code{\link{data.frame}} with gene information
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' result <- ensembleFS(data,
#'             decisions,
#'             methods = c('fs.utest', 'fs.mrmr'),
#'             method.cv = 'kfoldcv',
#'             params.cv = list(k = 3, iter = 10),
#'             level.cor = 0.75,
#'             params = list(adjust = 'fdr', feature.number = 10, alpha = 0.05),
#'             asm = c('fs.utest', 'fs.mrmr'),
#'             model = c('fs.utest', 'fs.mrmr')
#'             )
#'
#' gene.top <- get.top.gene(result$selected.feature, 15 , 20)
#' info.gene <- gene.info.top.gene(gene.top, condition.methods = 'union')
#'
#' }
#'
#' @export
get.info.top.gene <- function(gene.top, condition.methods = 'union'){
  if(condition.methods == 'intersect'){
    var.imp <- Reduce(intersect, gene.top)
  }
  else if(condition.methods == 'union'){
    var.imp <- unique(unlist(gene.top))
  }
  result <-  get.info.gprofiler(var.imp)
}
