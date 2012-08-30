\name{mvdcoord}
\alias{mvdcoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mean/variance differences discriminant coordinates}
\description{
  Discriminant projections as defined in Young, Marco and Odell (1987).
  The principle is to maximize the projection of a matrix consisting of
  the differences between the means of all classes and the first mean
  and the differences between the covariance matrices of all classes and
  the forst covariance matrix.
}
\usage{
mvdcoord(xd, clvecd, clnum=1, sphere="mcd", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{clnum}{integer. Number of the class to which all differences are
    computed.}
  \item{sphere}{a covariance matrix or one of
    "mve", "mcd", "classical", "none". The matrix used for sphering the
    data. "mcd" and "mve" are robust covariance matrices as implemented
    in \code{\link{cov.rob}}. "classical" refers to the classical
    covariance matrix. "none" means no sphering and use of the raw
    data.}
  \item{...}{no effect}
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
  Young, D. M., Marco, V. R. and Odell, P. L. (1987). Quadratic
  discrimination: some results on optimal low-dimensional
  representation, \emph{Journal of Statistical Planning and Inference},
  17, 307-319.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{plotcluster}} for straight forward discriminant plots.
  \code{\link{discrproj}} for alternatives.
  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
  set.seed(4634)
  face <- rFace(300,dMoNo=2,dNoEy=0,p=3)
  grface <- as.integer(attr(face,"grouping"))
  mcf <- mvdcoord(face,grface)
  plot(mcf$proj,col=grface)
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



