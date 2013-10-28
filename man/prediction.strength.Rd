\name{prediction.strength}
\alias{prediction.strength}
\alias{print.predstr}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction strength for estimating number of clusters}
\description{
  Computes the prediction strength of a clustering of a dataset into
  different numbers of components. The prediction strength is
  defined according to Tibshirani and Walther (2005), who recommend to
  choose as optimal number of cluster the largest number of clusters
  that leads to a prediction strength above 0.8 or 0.9. See details.

  Various clustering methods can be used, see argument \code{clustermethod}. 
}
\usage{
  prediction.strength(xdata, Gmin=2, Gmax=10, M=50,
                      clustermethod=kmeansCBI,
                                classification="centroid",
                                cutoff=0.8,nnk=1,
                      distances=inherits(xdata,"dist"),count=FALSE,...)
  \method{print}{predstr}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{data (something that can be coerced into a matrix). Note
    that this can currently not be a dissimilarity matrix.}
  \item{Gmin}{integer. Minimum number of clusters. Note that the
    prediction strength for 1 cluster is trivially 1, which is
    automatically included if \code{GMin>1}. Therefore \code{GMin<2} is
    useless.}  
  \item{Gmax}{integer. Maximum number of clusters.}
  \item{M}{integer. Number of times the dataset is divided into two
    halves.}
  \item{clustermethod}{an interface function (the function name, not a
    string containing the name, has to be provided!). This defines the
    clustering method. See the "Details"-section of \code{\link{clusterboot}}
    and \code{\link{kmeansCBI}} for the format. Clustering methods for
    \code{prediction.strength} must have a \code{k}-argument for the number of
    clusters, must operate on n times p data matrices
    and must otherwise follow the specifications in \code{\link{clusterboot}}.}
  \item{classification}{string.
    This determines how non-clustered points are classified to given
    clusters. Options are explained in \code{\link{classifnp}}.
    Certain classification methods are connected to certain clustering
    methods. \code{classification="averagedist"} is recommended for
    average linkage, \code{classification="centroid"} is recommended for
    k-means, clara and pam, \code{classification="knn"} with
    \code{nnk=1} is recommended for single linkage and
    \code{classification="qda"} is recommended for Gaussian mixtures
    with flexible covariance matrices.  
  }
  \item{cutoff}{numeric between 0 and 1. The optimal number of clusters
    is the maximum one with prediction strength above \code{cutoff}.}
  \item{nnk}{number of nearest neighbours if
    \code{classification="knn"}, see \code{\link{classifnp}}.}
  \item{distances}{logical. If \code{TRUE}, data will be interpreted as
    dissimilarity matrix, passed on to clustering methods as
    \code{"dist"}-object, and \code{\link{classifdist}} will be used for
    classification.}   
  \item{count}{logical. \code{TRUE} will print current number of
    clusters and simulation run number on the screen.}
  \item{x}{object of class \code{predstr}.}
  \item{...}{arguments to be passed on to the clustering method.}
}

\value{
  \code{prediction.strength} gives out an object of class
  \code{predstr}, which is a
  list with components
  \item{predcorr}{list of vectors of length \code{M} with relative
  frequencies of correct predictions (clusterwise minimum). Every list
  entry refers to a certain number of clusters.}
  \item{mean.pred}{means of \code{predcorr} for all numbers of
    clusters.}
  \item{optimalk}{optimal number of clusters.}
  \item{cutoff}{see above.}
  \item{method}{a string identifying the clustering method.}
  \item{Gmax}{see above.}
  \item{M}{see above.}
}

\details{
  The prediction strength for a certain number of clusters k under a
  random partition of the dataset in halves A and B is defined as
  follows. Both halves are clustered with k clusters. Then the points of
  A are classified to the clusters of B. In the original paper
  this is done by assigning every
  observation in A to the closest cluster centroid in B (corresponding
  to \code{classification="centroid"}), but other methods are possible,
  see \code{\link{classifnp}}. A pair of points A in
  the same A-cluster is defined to be correctly predicted if both points
  are classified into the same cluster on B. The same is done with the
  points of B relative to the clustering on A. The prediction strength
  for each of the clusterings is the minimum (taken over all clusters)
  relative frequency of correctly predicted pairs of points of that
  cluster. The final mean prediction strength statistic is the mean over
  all 2M clusterings.
}

\references{
  Tibshirani, R. and Walther, G. (2005) Cluster Validation by 
Prediction Strength, \emph{Journal of Computational and Graphical 
  Statistics}, 14, 511-528.
}

\seealso{
  \code{\link{kmeansCBI}}, \code{\link{classifnp}}
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  prediction.strength(iriss,2,3,M=3)
  prediction.strength(iriss,2,3,M=3,clustermethod=claraCBI)
# The examples are fast, but of course M should really be larger.
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



