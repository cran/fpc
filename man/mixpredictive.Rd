\name{mixpredictive}
\alias{mixpredictive}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction strength of merged Gaussian mixture}
\description{
  Computes the prediction strength of clustering by
  merging Gaussian mixture components, see \code{\link{mergenormals}}.
  The predictive strength is
  defined according to Tibshirani and Walther (2005), carried out as
  described in Hennig (2010), see details. 
}
\usage{
  mixpredictive(xdata, Gcomp, Gmix, M=50, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{data (something that can be coerced into a matrix).}
  \item{Gcomp}{integer. Number of components of the underlying Gaussian mixture.}
  \item{Gmix}{integer. Number of clusters after merging Gaussian components.}
  \item{M}{integer. Number of times the dataset is divided into two
    halves.} 
  \item{...}{further arguments that can potentially arrive in calls but
    are currently not used.}
}

\value{
  List with components
  \item{predcorr}{vector of length \code{M} with relative frequencies of
    correct predictions (clusterwise minimum).}
  \item{mean.pred}{mean of \code{predcorr}.}
}

\details{
  The prediction strength for a certain number of clusters \code{Gmix} under a
  random partition of the dataset in halves A and B is defined as
  follows. Both halves are clustered with \code{Gmix}
  clusters. Then the points of
  A are classified to the clusters of B. This is done by use of the
  maximum a posteriori rule for mixtures as in Hennig (2010),
  differently from Tibshirani and Walther (2005). A pair of points A in
  the same A-cluster is defined to be correctly predicted if both points
  are classified into the same cluster on B. The same is done with the
  points of B relative to the clustering on A. The prediction strength
  for each of the clusterings is the minimum (taken over all clusters)
  relative frequency of correctly predicted pairs of points of that
  cluster. The final mean prediction strength statistic is the mean over
  all 2M clusterings.
}

\references{
  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
  
  Tibshirani, R. and Walther, G. (2005) Cluster Validation by 
Prediction Strength, \emph{Journal of Computational and Graphical 
  Statistics}, 14, 511-528.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{prediction.strength}} for Tibshirani and Walther's
  original method.
  \code{\link{mergenormals}} for the clustering method applied here.
}

\examples{
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  mp <- mixpredictive(iriss,2,2,M=2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}

