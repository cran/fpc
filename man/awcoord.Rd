\name{awcoord}
\alias{awcoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Asymmetric weighted discriminant coordinates}
\description{
  Asymmetric weighted discriminant coordinates as defined
  in Hennig (2003). Asymmetric discriminant projection means that there
  are two classes, one of which is treated as the homogeneous class
  (i.e., it should appear homogeneous and separated in the resulting projection)
  while the other may be heterogeneous. 
  The principle is to maximize the ratio between the projection of a between
  classes separation matrix and the projection of the covariance matrix
  within the homogeneous class. Points are weighted according to their
  (robust) Mahalanobis distance to the homogeneous class. 
}
\usage{
awcoord(xd, clvecd, clnum=1, mahal="square", method="classical",
                     clweight=switch(method,classical=FALSE,TRUE),
                     alpha=0.99, subsample=0, countmode=1000, ...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{clnum}{integer. Number of the homogeneous class.}
  \item{mahal}{"md" or "square". If "md", the points are weighted by the
    square root of the \code{alpha}-quantile of the
    corresponding chi squared distribution
    over the roots of their Mahalanobis distance to the
    homogeneous class, unless
    this is smaller than 1. If "square" (which is recommended), the
    (originally squared) Mahalanobis distance and the
    unrooted quantile is used.}
  \item{method}{one of
    "mve", "mcd" or "classical". Covariance matrix used within the
    homogeneous class and for the computation of the Mahalanobis distances.
    "mcd" and "mve" are robust covariance matrices as implemented
    in \code{\link{cov.rob}}. "classical" refers to the classical
    covariance matrix.}
  \item{clweight}{logical. If \code{FALSE}, only the points of the
    heterogeneous class are weighted. This, together with
    \code{method="classical"}, computes AWC as defined in Hennig (2003). If
    \code{TRUE}, all points are weighted. This, together with
    \code{method="mcd"}, computes ARC as defined in Hennig (2003).}
  \item{alpha}{numeric between 0 and 1. The corresponding quantile of
    the chi squared distribution is used for the downweighting
    of points. Points with a smaller Mahalanobis distance to the
    homogeneous class get full weight.}
  \item{subsample}{integer. If 0, all points are used. Else, only a
    subsample of \code{subsample} of the points is used.}
  \item{countmode}{optional positive integer. Every \code{countmode}
    algorithm runs \code{awcoord} shows a message.}
  \item{...}{no effect}
}
\details{
  The square root of the homogeneous classes covariance matrix
  is inverted by use of
  \code{\link{tdecomp}}, which can be expected to give
  reasonable results for singular within-class covariance matrices.
}
% \details{
% }
\value{
  List with the following components
  \item{ev}{eigenvalues in descending order.}
  \item{units}{columns are coordinates of projection basis vectors.
  New points \code{x} can be projected onto the projection basis vectors
  by \code{x \%*\% units}}
  \item{proj}{projections of \code{xd} onto \code{units}.} 
}
\references{
  Hennig, C. (2003) Symmetric, asymmetric, and robust linear dimension
  reduction for classification, submitted,
  \url{http://stat.ethz.ch/Research-Reports/108.html}.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}
}

\seealso{
  \code{\link{plotcluster}} for straight forward discriminant plots.
  \code{\link{discrproj}} for alternatives.
  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
  set.seed(4634)
  face <- rFace(600,dMoNo=2,dNoEy=0)
  grface <- as.integer(attr(face,"grouping"))
  awcf <- awcoord(face,grface==1)
  # awcf2 <- ancoord(face,grface==1, method="mcd")
  plot(awcf$proj,col=1+(grface==1))
  # plot(awcf2$proj,col=1+(grface==1))
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



