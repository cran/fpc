\name{distrsimilarity}
\alias{distrsimilarity}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Similarity of within-cluster distributions to normal and uniform}
\description{
  Two measures of dissimilarity between the within-cluster distributions of
  a dataset and normal or uniform distribution. For the normal it's the
  Kolmogorov dissimilarity between the Mahalanobis distances to the
  center and a chi-squared distribution. For the uniform it is the
  Kolmogorov distance between the distance to the kth nearest neighbour
  and a Gamma distribution (this is based on Byers and Raftery (1998)).
  The clusterwise values are aggregated by weighting with the cluster sizes.
}
\usage{
distrsimilarity(x,clustering,noisecluster = FALSE,
distribution=c("normal","uniform"),nnk=2,
largeisgood=FALSE,messages=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clustering}{integer vector of class numbers; length must equal
    \code{nrow(x)}, numbers must go from 1 to the number of clusters.}
  \item{noisecluster}{logical. If \code{TRUE}, the cluster with the
    largest number is ignored for the computations.}
  \item{distribution}{vector of \code{"normal", "uniform"} or
    both. Indicates which of the two dissimilarities is/are computed.}
  \item{nnk}{integer. Number of nearest neighbors to use for
    dissimilarity to the uniform.}
  \item{largeisgood}{logical. If \code{TRUE}, dissimilarities are
    transformed to \code{1-d} (this means that larger values indicate a
    better fit).}
  \item{messages}{logical. If \code{TRUE}, warnings are given if
    within-cluster covariance matrices are not invertible (in which case
    all within-cluster Mahalanobis distances are set to zero).}
}
\note{
  It is very hard to capture similarity to a multivariate normal or
  uniform in a single value, and both used here have their
  shortcomings. Particularly, the dissimilarity to the uniform can still
  indicate a good fit if there are holes or it's a uniform distribution
  concentrated on several not connected sets.  
}
% \details{
% }
\value{
  List with the following components
  \item{kdnorm}{Kolmogorov distance between distribution of
    within-cluster Mahalanobis
    distances and appropriate chi-squared distribution, aggregated over
    clusters (I am grateful to Agustin Mayo-Iscar for the idea).}
  \item{kdunif}{Kolmogorov distance between distribution of distances to
    \code{nnk}th nearest within-cluster neighbor and appropriate
    Gamma-distribution, see Byers and Raftery (1998), aggregated over
    clusters.}
  \item{kdnormc}{vector of cluster-wise Kolmogorov distances between
    distribution of within-cluster Mahalanobis
    distances and appropriate chi-squared distribution.}
  \item{kdunifc}{vector of cluster-wise Kolmogorov distances between
    distribution of distances to \code{nnk}th nearest within-cluster
    neighbor and appropriate Gamma-distribution.}
  \item{xmahal}{vector of Mahalanobs distances to the respective cluster
    center.}
  \item{xdknn}{vector of distance to \code{nnk}th nearest within-cluster
    neighbor.}
}
\references{
Byers, S. and Raftery, A. E. (1998) Nearest-Neighbor Clutter
Removal for Estimating Features in Spatial Point Processes,
\emph{Journal of the American Statistical Association}, 93, 577-584.

Hennig, C. (2017) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Proceedings of ASMDA 2017}, 501-520,
    \url{https://arxiv.org/abs/1703.09282}


}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{cqcluster.stats}},\code{\link{cluster.stats}}
  for more cluster validity statistics.
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  km3 <- kmeans(face,3)
  distrsimilarity(face,km3$cluster) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line
\keyword{cluster}% __ONLY ONE__ keyword per line



