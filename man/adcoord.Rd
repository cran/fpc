\name{adcoord}
\alias{adcoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Asymmetric discriminant coordinates}
\description{
  Asymmetric discriminant coordinates as defined
  in Hennig (2003). Asymmetric discriminant projection means that there
  are two classes, one of which is treated as the homogeneous class
  (i.e., it should appear homogeneous and separated in the resulting projection)
  while the other may be heterogeneous. 
  The principle is to maximize the ratio between the projection of a between
  classes separation matrix and the projection of the covariance matrix
  within the homogeneous class.  
}
\usage{
adcoord(xd, clvecd, clnum=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{clnum}{integer. Number of the homogeneous class.}
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
  adcf <- adcoord(face,grface==2)
  adcf2 <- adcoord(face,grface==4)
  plot(adcf$proj,col=1+(grface==2))
  plot(adcf2$proj,col=1+(grface==4))
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



