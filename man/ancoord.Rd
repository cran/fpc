\name{ancoord}
\alias{ancoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Asymmetric neighborhood based discriminant coordinates}
\description{
  Asymmetric neighborhood based discriminant coordinates as defined
  in Hennig (2003). Asymmetric discriminant projection means that there
  are two classes, one of which is treated as the homogeneous class
  (i.e., it should appear homogeneous and separated in the resulting projection)
  while the other may be heterogeneous. 
  The principle is to maximize the ratio between the projection of a between
  classes covariance matrix, which is defined by averaging the
  between classes covariance matrices in the neighborhoods of the points
  of the homogeneous class and the projection of the covariance matrix
  within the homogeneous class.  
}
\usage{
ancoord(xd, clvecd, clnum=1, nn=50, method="mcd", countmode=1000, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{clnum}{integer. Number of the homogeneous class.}
  \item{nn}{integer. Number of points which belong to the neighborhood
    of each point (including the point itself).}
  \item{method}{one of
    "mve", "mcd" or "classical". Covariance matrix used within the
    homogeneous class.
    "mcd" and "mve" are robust covariance matrices as implemented
    in \code{\link{cov.rob}}. "classical" refers to the classical
    covariance matrix.}
  \item{countmode}{optional positive integer. Every \code{countmode}
    algorithm runs \code{ancoord} shows a message.}
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
  Hennig, C. (2004) Asymmetric linear dimension reduction for classification.
  Journal of Computational and Graphical Statistics 13, 930-945 .
  
  Hennig, C. (2005)  A method for visual cluster validation.  In:
  Weihs, C. and Gaul, W. (eds.): Classification - The Ubiquitous
  Challenge. Springer, Heidelberg 2005, 153-160.

}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
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
  ancf2 <- ancoord(face,grface==4)
  plot(ancf2$proj,col=1+(grface==4))
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



