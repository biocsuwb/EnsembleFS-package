% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs.mcfs.R
\name{fs.mcfs}
\alias{fs.mcfs}
\title{Monte Carlo Feature Selection}
\usage{
fs.mcfs(x, y, params = list(cutoff.method = c("kmeans")))
}
\arguments{
\item{x}{input data where columns are variables and rows are observations (all numeric)}

\item{y}{decision variable as a boolean vector of length equal to number of observations}

\item{params}{A \code{\link{list}} with the following fields:
\itemize{
\item \code{cutoff.method} -- cutoff method MCFS: "permutations", "criticalAngle", "kmeans", "mean", "contrast"
}}
}
\value{
A \code{\link{data.frame}} with selected features
}
\description{
Monte Carlo Feature Selection
}
\details{
Performs Monte Carlo Feature Selection (MCFS) on a given data set.
}
\examples{
\dontrun{

decisions <- data$class
data$class <- NULL

fs.mcfs(data, decisions, ...)
}

}
