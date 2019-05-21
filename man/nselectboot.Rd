\name{nselectboot}
\alias{nselectboot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Selection of the number of clusters via bootstrap}
\description{
  Selection of the number of clusters via bootstrap as explained in Fang
  and Wang (2012). Several times 2 bootstrap samples are drawn from the
  data and the number of clusters is chosen by optimising an instability
  estimation from these pairs.

  In principle all clustering methods can be used that have a
  CBI-wrapper, see \code{\link{clusterboot}},
  \code{\link{kmeansCBI}}. However, the currently implemented
  classification methods are not necessarily suitable for all of them,
  see argument \code{classification}.
}
\usage{
nselectboot(data,B=50,distances=inherits(data,"dist"),
                        clustermethod=NULL,
                        classification="averagedist",krange=2:10,
                        count=FALSE,nnk=1, ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{something that can be coerced into a matrix. The data
    matrix - either an \code{n*p}-data matrix (or data frame) or an
    \code{n*n}-dissimilarity matrix (or \code{dist}-object).}
  \item{B}{integer. Number of resampling runs.}
  \item{distances}{logical. If \code{TRUE}, the data is interpreted as
    dissimilarity matrix. If \code{data} is a \code{dist}-object,
    \code{distances=TRUE} automatically, otherwise
    \code{distances=FALSE} by default. This means that you have to set
    it to \code{TRUE} manually if \code{data} is a dissimilarity matrix.}
  \item{clustermethod}{an interface function (the function name, not a
    string containing the name, has to be provided!). This defines the
    clustering method. See the "Details"-section of \code{\link{clusterboot}}
    and \code{\link{kmeansCBI}} for the format. Clustering methods for
    \code{nselectboot} must have a \code{k}-argument for the number of
    clusters and must otherwise follow the specifications in
    \code{\link{clusterboot}}. Note that \code{nselectboot} won't work
    with CBI-functions that implicitly already estimate the number of
    clusters such as \code{\link{pamkCBI}}; use \code{\link{claraCBI}}
    if you want to run it for pam/clara clustering. 
  }
  \item{classification}{string.
    This determines how non-clustered points are classified to given
    clusters. Options are explained in \code{\link{classifdist}} (if
    \code{distances=TRUE}) and \code{\link{classifnp}} (otherwise).
    Certain classification methods are connected to certain clustering
    methods. \code{classification="averagedist"} is recommended for
    average linkage, \code{classification="centroid"} is recommended for
    k-means, clara and pam (with distances it will work with
    \code{\link{claraCBI}} only), \code{classification="knn"} with
    \code{nnk=1} is recommended for single linkage and
    \code{classification="qda"} is recommended for Gaussian mixtures
    with flexible covariance matrices.  
  }
  \item{krange}{integer vector; numbers of clusters to be tried.}
  \item{count}{logical. If \code{TRUE}, numbers of clusters and
    bootstrap runs are printed.}
  \item{nnk}{number of nearest neighbours if
    \code{classification="knn"}, see  \code{\link{classifdist}} (if
    \code{distances=TRUE}) and \code{\link{classifnp}} (otherwise).}
  \item{...}{arguments to be passed on to the clustering method.}
}
    
\value{
  \code{nselectboot} returns a list with components
  \code{kopt,stabk,stab}.
  \item{kopt}{optimal number of clusters.}
  \item{stabk}{mean instability values for numbers of clusters.}
  \item{stab}{matrix of instability values for all bootstrap runs and
    numbers of clusters.} 
}
\references{
  Fang, Y. and Wang, J. (2012) Selection of the number of clusters via
  the bootstrap method. \emph{Computational Statistics and Data
    Analysis}, 56, 468-477.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{classifdist}}, \code{\link{classifnp}},
  \code{\link{clusterboot}},\code{\link{kmeansCBI}}
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  nselectboot(dist(face),B=2,clustermethod=disthclustCBI,
   method="average",krange=5:7)
  nselectboot(dist(face),B=2,clustermethod=claraCBI,
   classification="centroid",krange=5:7)
  nselectboot(face,B=2,clustermethod=kmeansCBI,
   classification="centroid",krange=5:7)
# Of course use larger B in a real application.
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



