\name{dbscan}
\alias{dbscan}
\alias{print.dbscan}
\alias{plot.dbscan}
\alias{predict.dbscan}
\title{DBSCAN density reachability and connectivity clustering}
\description{
  Generates a density based clustering of arbitrary shape as introduced
  in Ester et al. (1996).
}
\usage{
  dbscan(data, eps, MinPts = 5, scale = FALSE, method = c("hybrid", "raw",
    "dist"), seeds = TRUE, showplot = FALSE, countmode = NULL)
  \method{print}{dbscan}(x, ...)
  \method{plot}{dbscan}(x, data, ...)
  \method{predict}{dbscan}(object, data, newdata = NULL,
predict.max=1000, ...)
}
\arguments{
  \item{data}{data matrix, data.frame, dissimilarity matrix or
    \code{dist}-object. Specify \code{method="dist"} if the data should
    be interpreted as dissimilarity matrix or object. Otherwise
    Euclidean distances will be used.}
  \item{eps}{ Reachability distance, see Ester et al. (1996). }
  \item{MinPts}{ Reachability minimum no. of points, see Ester et al. (1996). }
  \item{scale}{ scale the data if \code{TRUE}. }
  \item{method}{ "dist" treats data as distance matrix (relatively fast
    but memory expensive), "raw" treats data as raw data and avoids
    calculating a distance matrix (saves memory but may be slow),
    "hybrid" expects also raw data, but calculates partial distance
    matrices (very fast with moderate memory requirements).}
  \item{seeds}{FALSE to not include the \code{isseed}-vector in the
    \code{dbscan}-object.}
  \item{showplot}{ 0 = no plot, 1 = plot per iteration, 2 = plot per
    subiteration. }
  \item{countmode}{ NULL or vector of point numbers at which to report
    progress. }
  \item{x}{object of class \code{dbscan}.}
  \item{object}{object of class \code{dbscan}.}
  \item{newdata}{ matrix or data.frame with raw data to predict. }
  \item{predict.max}{ max. batch size for predictions. } 
  \item{...}{Further arguments transferred to plot methods.}
}
\details{
  Clusters require a minimum no of points (MinPts) within a maximum distance
  (eps) around one of its members (the seed).
  Any point within eps around any point which satisfies the seed condition
  is a cluster member (recursively).
  Some points may not belong to any clusters (noise).
  
  We have clustered a 100.000 x 2 dataset in 40 minutes on a Pentium M 1600
  MHz.

  \code{print.dbscan} shows a statistic of the number of points
  belonging to the clusters that are seeds and border points.

  \code{plot.dbscan} distinguishes between seed and border points by
  plot symbol.

}
\value{
  \code{predict.dbscan} gives out a vector of predicted clusters for the
  points in \code{newdata}.
  
  \code{dbscan} gives out 
  an object of class 'dbscan' which is a LIST with components
  \item{cluster}{integer vector coding cluster membership with noise
    observations (singletons) coded as 0 }
  \item{isseed}{logical vector indicating whether a point is a seed (not
    border, not noise)}
  \item{eps}{parameter eps}
  \item{MinPts}{parameter MinPts}
}
\references{ Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu
(1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial
Databases with Noise. Institute for Computer Science, University of Munich.
Proceedings of 2nd International Conference on Knowledge Discovery and Data
Mining (KDD-96). }
\author{Jens Oehlschlaegel, based on a draft by Christian Hennig.}
\note{ this is a simplified version of the original algorithm (no K-D-trees
  used), thus we have \eqn{o(n^2)} instead of \eqn{o(n*log(n))} }
\examples{
  set.seed(665544)
  n <- 600
  x <- cbind(runif(10, 0, 10)+rnorm(n, sd=0.2), runif(10, 0, 10)+rnorm(n,
    sd=0.2))
  par(bg="grey40")
  ds <- dbscan(x, 0.2)
# run with showplot=1 to see how dbscan works.
  ds
  plot(ds, x)

  x2 <- matrix(0,nrow=4,ncol=2)
  x2[1,] <- c(5,2)
  x2[2,] <- c(8,3)
  x2[3,] <- c(4,4)
  x2[4,] <- c(9,9)
  predict(ds, x, x2)

  n <- 600
  x <- cbind((1:3)+rnorm(n, sd=0.2), (1:3)+rnorm(n, sd=0.2))

# Not run, but results from my machine are 0.105 - 0.068 - 0.255:
#  system.time(ds <- dbscan(x, 0.3, countmode=NULL, method="raw"))[3] 
#  system.time(dsb <- dbscan(x, 0.3, countmode=NULL, method="hybrid"))[3]
#  system.time(dsc <- dbscan(dist(x), 0.3, countmode=NULL,
#    method="dist"))[3]
}
\keyword{multivariate}
\keyword{cluster}

