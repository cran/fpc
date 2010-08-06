\name{kmeansruns}
\alias{kmeansruns}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{k-means with estimating k and initialisations}
\description{
  This calls the function \code{\link{kmeans}} to perform a k-means
  clustering, but initializes the k-means algorithm several times with
  random points from the data set as means. Furthermore, it is more
  robust against the occurrence of empty clusters in the algorithm and
  it estimates the number of clusters by either the Calinski Harabasz
  index (\code{\link{calinhara}}) or average silhouette width (see
  \code{\link[cluster]{pam.object}}). The Duda-Hart test
  (\code{\link{dudahart2}}) is applied to decide whether there should be
  more than one cluster (unless 1 is excluded as number of clusters).
}
\usage{
kmeansruns(data,krange=2:10,criterion="ch",
                       iter.max=100,runs=100,
                       scaledata=FALSE,alpha=0.001,
                       critout=FALSE,plot=FALSE,...)
}
\arguments{
  \item{data}{A numeric matrix of data, or an object that can be coerced to
          such a matrix (such as a numeric vector or a data frame with
          all numeric columns). }
  \item{krange}{integer vector. Numbers of clusters which are to be
    compared by the average silhouette width criterion. Note: average
    silhouette width and Calinski-Harabasz can't estimate number of
    clusters \code{nc=1}. If 1 is included, a Duda-Hart test is applied
    and 1 is estimated if this is not significant.}
  \item{criterion}{one of \code{"asw"} or \code{"ch"}. Determines
    whether average silhouette width or Calinski-Harabasz is applied.}
  \item{iter.max}{integer. The maximum number of iterations allowed.}
  \item{runs}{integer. Number of starts of the k-means algorithm.}
  \item{scaledata}{logical. If \code{TRUE}, the variables are centered
    and scaled to unit variance before execution.}
  \item{alpha}{numeric between 0 and 1, tuning constant for
    \code{\link{dudahart2}} (only used for 1-cluster test).}
  \item{critout}{logical. If \code{TRUE}, the criterion value is printed
    out for every number of clusters.}
  \item{plot}{logical. If \code{TRUE}, every clustering resulting from a
    run of the algorithm is plotted.}
  \item{...}{further arguments to be passed on to \code{\link{kmeans}}.}
}

\value{
  The output of the optimal run of the \code{\link{kmeans}}-function.
  A list with components
  \item{cluster}{A vector of integers indicating the cluster to which each
          point is allocated.} 
  \item{centers}{A matrix of cluster centers.}
  \item{withinss}{The within-cluster sum of squares for each cluster.}
  \item{size}{The number of points in each cluster.}
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\references{
    Calinski, R. B., and Harabasz, J. (1974) A Dendrite Method for Cluster 
Analysis, \emph{Communications in Statistics}, 3, 1-27.

  Duda, R. O. and Hart, P. E. (1973) \emph{Pattern Classification and
  Scene Analysis}. Wiley, New York.

Hartigan, J. A. and Wong, M. A. (1979).  A K-means clustering
algorithm. \emph{Applied Statistics}, 28, 100-108.

     Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
     An Introduction to Cluster Analysis". Wiley, New York.
}
\seealso{
  \code{\link{kmeans}}, \code{\link{pamk}},
  \code{\link{calinhara}}, \code{\link{dudahart2}})  
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  pka <- kmeansruns(face,krange=1:5,critout=TRUE,runs=2,criterion="asw")
  pkc <- kmeansruns(face,krange=1:5,critout=TRUE,runs=2,criterion="ch")
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



