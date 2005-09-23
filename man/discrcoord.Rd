\name{discrcoord}
\alias{discrcoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Discriminant coordinates/canonical variates}
\description{
  Computes discriminant coordinates, sometimes referred to as "canonical
  variates" as described in Seber (1984).
}
\usage{
discrcoord(xd, clvecd, pool = "n", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{pool}{string. Determines how the within classes
    covariance is pooled. "n" means that the class covariances are
    weighted corresponding to the number of points in each class
    (default). "equal" means that all classes get equal weight.}
  \item{...}{no effect}
}
\details{
  The matrix T (see Seber (1984), p. 270) is inverted by use of
  \code{\link{tdecomp}}, which can be expected to give
  reasonable results for singular within-class covariance matrices.
}
\value{
  List with the following components
  \item{ev}{eigenvalues in descending order.}
  \item{units}{columns are coordinates of projection basis vectors.
  New points \code{x} can be projected onto the projection basis vectors
  by \code{x \%*\% units}}
  \item{proj}{projections of \code{xd} onto \code{units}.} 
}
\references{
  Seber, G. A. F. (1984). \emph{Multivariate Observations}. New York: Wiley.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{plotcluster}} for straight forward discriminant plots.
  
  \code{\link{batcoord}} for discriminating projections for two classes,
  so that also the differences in variance are shown (\code{discrcoord} is
  based only on differences in mean).

  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
  set.seed(4634)
  face <- rFace(600,dMoNo=2,dNoEy=0)
  grface <- as.integer(attr(face,"grouping"))
  dcf <- discrcoord(face,grface)
  plot(dcf$proj,col=grface)
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



