kfoldcv <- function(x, y, k){
  index.train <- list()
  index.test <- list()
  x$class <- y
  rownames(x) <- 1:nrow(x)
  as.numeric(row.names(x))
  data.Class0 <- as.data.frame(x[x$class == 0, 1:3])
  data.Class1 <- as.data.frame(x[x$class == 1, 1:3])
  size.Class0 <- nrow(data.Class0)
  size.Class1 <- nrow(data.Class1)
  index.permut.class0 <- sample.int(size.Class0 , size = size.Class0 , replace = FALSE)
  index.permut.class1 <- sample.int(size.Class1 , size = size.Class1 , replace = FALSE)
  data.Class0 <- data.Class0[index.permut.class0,]
  data.Class1 <- data.Class1[index.permut.class1,]
  sublist.kfold0 <- cut(seq(1,length(data.Class0[,1])), breaks = k , labels = FALSE)
  sublist.kfold1 <- cut(seq(1,length(data.Class1[,1])), breaks = k , labels = FALSE)
  for(k in 1:k){
    test.index.class0 <- which(sublist.kfold0 ==  k , arr.ind = TRUE)
    test.index.class1 <- which(sublist.kfold1 ==  k , arr.ind = TRUE)
    index.test[[k]] <- as.numeric(row.names(rbind(data.Class0[test.index.class0,] , data.Class1[test.index.class1,])))
    index.train[[k]] <- as.numeric(row.names(rbind(data.Class0[-test.index.class0,],data.Class1[-test.index.class1,])))

  }
  validation.index <- list(index.test, index.train)
  names(validation.index) <- c('testing', 'training')
  return(validation.index)
}


rsampling <- function(x, y, test.size = 0.3){
  index.train <- list()
  index.test <- list()
  size <- nrow(x)
  test.size <- as.integer(size * test.size)
  test.idx <- c()
  indexes <- 1:size
  for(i in 1:test.size){
    test.idx <- append(test.idx, sample(indexes, size = 1))
  }
  test.idx <- unique(test.idx)
  index.train[[1]] <- indexes[-test.idx]
  index.test[[1]] <- test.idx
  validation.index <- list(index.test, index.train)
  names(validation.index) <- c('testing', 'training')
  return(validation.index)
}



#' Cross-validation
#'
#' @details
#' creates a list with division of observational indices of the dataframe of
#' observational indices into training and test subsamples of a given number of iterations
#' Implements two methods : cross-validation k-fold and random sampling
#'
#' @param x input data where columns are variables and rows are observations (all numeric)
#' @param y decision variable as a boolean vector of length equal to number of observations
#' @param method cross-validation method \code{kfoldcv} for cross-validation \code{k-fold} or \code{rsampling} for \code{random sampling}
#' @param params.cv A \code{\link{list}} with the following fields:
#' \itemize{
#'   \item \code{niter} -- the number of validation repetitions
#'   \item \code{k} -- the number of groups that a given data sample is to be split into
#'   \item \code{test.size} -- testing set size for random sampling validation
#' }
#' @return A \code{\link{list}} with lindices divided into train and test subsamples niter repetitions.
#'
#' @examples
#' \dontrun{
#'
#' decisions <- data$class
#' data$class <- NULL
#'
#' list.index.cross <- cross.val(x = data,
#'                              y = decisions,
#'                              method = 'kfoldcv',
#'                              params.cv = list(niter = 10, k = 3))
#' }
#'
#' @export
cross.val <- function(x, y, method, params.cv = list(niter = 10, k = 3, test.size = 0.3)){
  if(!(method %in% c('kfoldcv', 'rsampling'))){
    stop('Unknown validation method, use please kfoldcv or rsampling')
  }
  indexes.cross.val <- list()
  if(method == 'kfoldcv'){
    if(params.cv$k < 2){
      stop('k < 2')
    }
    for(i in 1:params.cv$niter){
      index <-kfoldcv(x, y, params.cv$k)
      indexes.cross.val <- append(indexes.cross.val, list(index))
    }
  }
  else if(method == 'rsampling'){
    if(params.cv$test.size > 0.9 || params.cv$test.size < 0.05){stop('Invalid test size for Random Sampling')}
    for(i in 1:params.cv$niter){
      index <-rsampling(x, y, params.cv$test.size)
      indexes.cross.val <- append(indexes.cross.val, list(index))
    }
  }
  return(indexes.cross.val)
}



