\name{solvecov}
\alias{solvecov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Inversion of (possibly singular) symmetric matrices}
\description{
  Tries to invert a matrix by \code{solve}. If this fails because of
  singularity, an
  eigenvector decomposition is computed, and eigenvalues below
  \code{1/cmax} are replaced by \code{1/cmax}, i.e., \code{cmax} will be
  the corresponding eigenvalue of the inverted matrix.
}
\usage{
solvecov(m, cmax = 1e+10)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{m}{a numeric symmetric matrix.}
  \item{cmax}{a positive value, see above.}
}

\value{
  A list with the following components:
  \item{inv}{the inverted matrix}
  \item{coll}{\code{TRUE} if \code{solve} failed because of singularity.}
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{\code{\link{solve}}, \code{\link{eigen}}}

\examples{
  x <- c(1,0,0,1,0,1,0,0,1)
  dim(x) <- c(3,3)
  solvecov(x)
}
\keyword{array}% at least one, from doc/KEYWORDS


