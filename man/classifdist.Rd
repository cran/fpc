\name{classifdist}
\alias{classifdist}
\alias{classifnp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Classification of unclustered points}
\description{
  Various methods for classification of unclustered points from
  clustered points for use within functions \code{nselectboot}
  and \code{prediction.strength}. 
}
\usage{
classifdist(cdist,clustering,
                      method="averagedist",
                      centroids=NULL,nnk=1)

classifnp(data,clustering,
                      method="centroid",cdist=NULL,
                      centroids=NULL,nnk=1)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{cdist}{dissimilarity matrix or \code{dist}-object. Necessary for
    \code{classifdist} but optional for \code{classifnp} and there only
    used if \code{method="averagedist"} (if not provided, \code{dist} is
    applied to \code{data}).}
  \item{data}{something that can be coerced into a an
    \code{n*p}-data matrix.}
  \item{clustering}{integer vector. Gives the cluster number (between 1
    and k for k clusters) for clustered points and should be -1 for
    points to be classified.}
  \item{method}{one of \code{"averagedist", "centroid", "qda",
      "knn"}. See details.}
  \item{centroids}{for \code{classifnp} a k times p matrix of cluster
    centroids. For \code{classifdist} a vector of numbers of centroid
    objects as provided by \code{\link[cluster]{pam}}. Only used if
    \code{method="centroid"}; in that case mandatory for
    \code{classifdist} but optional for \code{classifnp}, where cluster mean
    vectors are computed if \code{centroids=NULL}.}
  \item{nnk}{number of nearest neighbours if \code{method="knn"}.}
}

\details{
  \code{classifdist} is for data given as dissimilarity matrix,
  \code{classifnp} is for data given as n times p data matrix. 
  The following methods are supported:
  \describe{
    \item{"centroid"}{assigns observations to the cluster with closest
      cluster centroid as specified in argument \code{centroids} (this
      is associated to k-means and pam/clara-clustering).}
    \item{"averagedist"}{assigns to the cluster to which an observation
      has the minimum average dissimilarity to all points in the cluster
      (this is associated to average linkage clustering).}
    \item{"qda"}{only in \code{classifnp}. Classifies by quadratic
      discriminant analysis (this is associated to Gaussian clusters
      with flexible covariance matrices), calling
      \code{\link[MASS]{qda}} with default settings. If
      \code{\link[MASS]{qda}} gives an error (usually because a class
      was too small), \code{\link[MASS]{lda}} is used.}
    \item{"knn"}{classifies by \code{nnk} nearest neighbours (for
      \code{nnk=1}, this is associated with single linkage clustering).
      Calls \code{\link[class]{knn}} in \code{classifnp}.}
  }
}
    
\value{
  An integer vector giving cluster numbers for all observations; those
  for the observations already clustered in the input are the same as in
  the input.
}
  
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{prediction.strength}}, \code{\link{nselectboot}}
}
\examples{  
set.seed(20000)
x1 <- rnorm(50)
y <- rnorm(100)
x2 <- rnorm(40,mean=20)
x3 <- rnorm(10,mean=25,sd=100)
x <-cbind(c(x1,x2,x3),y)
truec <- c(rep(1,50),rep(2,40),rep(3,10))
topredict <- c(1,2,51,52,91)
clumin <- truec
clumin[topredict] <- -1

classifnp(x,clumin, method="averagedist")
classifnp(x,clumin, method="qda")
classifdist(dist(x),clumin, centroids=c(3,53,93),method="centroid")
classifdist(dist(x),clumin,method="knn")

}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



