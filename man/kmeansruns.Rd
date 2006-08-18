\name{kmeansruns}
\alias{kmeansruns}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{k-means clustering with several random initializations}
\description{
  This calls the function \code{\link{kmeans}} to perform a k-means
  clustering, but initializes the k-means algorithm several times with
  random points from the data set as means. Furthermore, it is more
  robust against the occurrence of empty clusters in the algorithm.
}
\usage{
kmeansruns(data,k,iter.max=100,runs=100,scaledata=FALSE,plot=FALSE)
}
\arguments{
  \item{data}{A numeric matrix of data, or an object that can be coerced to
          such a matrix (such as a numeric vector or a data frame with
          all numeric columns). }
  \item{k}{integer. The number of clusters.}
  \item{iter.max}{integer. The maximum number of iterations allowed.}
  \item{runs}{integer. Number of starts of the k-means algorithm.}
  \item{scaledata}{logical. If \code{TRUE}, the variables are centered
    and scaled to unit variance before execution.}
  \item{plot}{logical. If \code{TRUE}, every clustering resulting from a
    run of the algorithm is plotted.} 
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
\seealso{
  \code{\link{kmeans}}
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  kmr1 <- kmeansruns(face,k=5,runs=1)
  kmr5 <- kmeansruns(face,k=5,runs=5)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



