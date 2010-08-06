\name{rFace}
\alias{rFace}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{"Face-shaped" clustered benchmark datasets}
\description{
  Generates "face-shaped" clustered benchmark datasets.
}
\usage{
rFace(n, p = 6, nrep.top = 2, smile.coef = 0.6, dMoNo = 1.2, dNoEy = 1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{integer greater or equal to 10. Number of points.}
  \item{p}{integer greater or equal to 2. Dimension.}
  \item{nrep.top}{integer. Number of repetitions of the hair-top point.}
  \item{smile.coef}{numeric. Coefficient for quadratic term used for
    generation of mouth-points. Positive values=>smile.}
  \item{dMoNo}{number. Distance from mouth to nose.}
  \item{dNoEy}{number. Minimum vertical distance from mouth to eyes.}
}
\details{
  The function generates a nice benchmark example for cluster
  analysis. 
  There are six "clusters" in this data, of which the first five are
  clearly homogeneous patterns, but with different distributional
  shapes and different qualities of separation. The clusters are
  distinguished only in the first two dimensions. The attribute
  \code{grouping} is a factor giving the cluster numbers, see below.
  The sixth group of
  points corresponds to some hairs, and is rather a collection of
  outliers than a cluster in itself. This group contains
  \code{nrep.top+2} points. Of the remaining points, 20\% belong to
  cluster 1, the chin (quadratic function plus noise).
  10\% belong to cluster 2, the right eye (Gaussian). 30\% belong to
  cluster 3, the mouth (Gaussian/squared Gaussian). 
  20\% belong to cluster 4, the nose (Gaussian/gamma), and
  20\% belong to cluster 5, the left eye (uniform).

  The distributions of the further
  variables are homogeneous over
  all points. The third dimension is exponentially distributed, the
  fourth dimension is Cauchy distributed, all further distributions are
  Gaussian.
  
  Please consider the source code for exact generation of the clusters.
}
\value{
  An \code{n} times \code{p} numeric matrix with attributes
  \item{grouping}{a factor giving the cluster memberships of the points.}
  \item{indexlist}{a list of six vectors containing the indices of points
    belonging to the six groups.}
}

\author{Martin Maechler
  \email{maechler@stat.math.ethz.ch}
  \url{http://stat.ethz.ch/~maechler/}
  
  Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\examples{
  set.seed(4634)
  face <- rFace(600,dMoNo=2,dNoEy=0)
  grface <- as.integer(attr(face,"grouping"))
  plot(face, col = grface)
#  pairs(face, col = grface, main ="rFace(600,dMoNo=2,dNoEy=0)")
}
\keyword{data}% at least one, from doc/KEYWORDS
