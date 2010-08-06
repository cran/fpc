\name{prediction.strength}
\alias{prediction.strength}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction strength for estimating number of clusters}
\description{
  Computes the prediction strength of a clustering of a dataset into
  different numbers of components. The predictive strength is
  defined according to Tibshirani and Walther (2005), who recommend to
  choose as optimal number of cluster the largest number of clusters
  that leads to a prediction strength above 0.8 or 0.9. See details.
}
\usage{
  prediction.strength(xdata, Gmin=2, Gmax=10, method="kmeans", M=50,
                                cutoff=0.8,...)
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
  \item{method}{one of \code{"kmeans"}, \code{"pam"} or \code{"clara"},
    specifying the clustering method to be applied.}
  \item{M}{integer. Number of times the dataset is divided into two
    halves.}
  \item{cutoff}{numeric between 0 and 1. The optimal number of clusters
    is the maximum one with prediction strength above \code{cutoff}.} 
  \item{...}{arguments to be passed on to the clustering method.}
}

\value{
  List with components
  \item{predcorr}{list of vectors of length \code{M} with relative
  frequencies of correct predictions (clusterwise minimum). Every list
  entry refers to a certain number of clusters.}
  \item{mean.pred}{means of \code{predcorr} for all numbers of
    clusters.}
  \item{optimalk}{optimal number of clusters.}
  \item{cutoff}{see above.}
}

\details{
  The prediction strength for a certain number of clusters k under a
  random partition of the dataset in halves A and B is defined as
  follows. Both halves are clustered with k clusters. Then the points of
  A are classified to the clusters of B. This is done by assigning every
  observation in A to the closest cluster centroid in B (using the
  function \code{\link[class]{knn1}}). A pair of points A in
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
  \code{\link{kmeans}}, \code{\link[cluster]{pam}},
  \code{\link[cluster]{clara}}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  prediction.strength(iriss,2,3,M=3)
  prediction.strength(iriss,2,3,M=3,method="pam")
# The examples are fast, but of course M should really be larger.
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



