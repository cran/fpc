\name{mahalanofix}
\alias{mahalanofix}
\alias{mahalanofuz}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mahalanobis distances from center of indexed points}
\description{
  Computes the vector of (classical or robust)
  Mahalanobis distances of all points of \code{x}
  to the center of the points indexed (or weighted)
  by \code{gv}. The latter also determine
  the covariance matrix.

  Thought for use within \code{\link{fixmahal}}.
}
\usage{
mahalanofix(x, n = nrow(as.matrix(x)), p = ncol(as.matrix(x)), gv =
rep(1, times = n), cmax = 1e+10, method = "ml")

mahalanofuz(x, n = nrow(as.matrix(x)), p = ncol(as.matrix(x)),
                         gv = rep(1, times=n), cmax = 1e+10) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a numerical data matrix, rows are points, columns are variables.}
  \item{n}{positive integer. Number of points.}
  \item{p}{positive integer. Number of variables.}
  \item{gv}{for \code{mahalanofix}
    a logical or 0-1 vector of length \code{n}. For \code{mahalanofuz} a
    numerical vector with values between 0 and 1.}
  \item{cmax}{positive number. used in \code{\link{solvecov}} if
    covariance matrix is singular.}
  \item{method}{\code{"ml"}, \code{"classical"},
    \code{"mcd"} or \code{"mve"}. Method to compute the covariance
    matrix estimator. See \code{\link{cov.rob}}, \code{\link{fixmahal}}.}
  
}
\details{
  \code{\link{solvecov}} is used to invert the covariance matrix. The methods
  \code{"mcd"} and \code{"mve"} in \code{mahalanofix} do not work properly
  with point constellations with singular covariance matrices!
}
\value{
  A list of the following components:
  \item{md}{vector of Mahalanobis distances.}
  \item{mg}{mean of the points indexed by \code{gv}, weighted mean in
    \code{mahalanofuz}.}
  \item{covg}{covariance matrix of the points indexed by \code{gv},
    weighted covariance matrix in \code{mahalanofuz}.} 
  \item{covinv}{\code{covg} inverted by \code{\link{solvecov}}.}
  \item{coll}{logical. If \code{TRUE}, \code{covg} has been
    (numerically) singular.}
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\note{Methods \code{"mcd"} and \code{"mve"} require library \code{lqs}.}

\seealso{\code{\link{fixmahal}}, \code{\link{solvecov}}, \code{\link{cov.rob}}}

\examples{
  x <- c(1,2,3,4,5,6,7,8,9,10)
  y <- c(1,2,3,8,7,6,5,8,9,10)
  mahalanofix(cbind(x,y),gv=c(0,0,0,1,1,1,1,1,0,0))
  mahalanofix(cbind(x,y),gv=c(0,0,0,1,1,1,1,0,0,0))
  mahalanofix(cbind(x,y),gv=c(0,0,0,1,1,1,1,1,0,0),method="mcd")
  mahalanofuz(cbind(x,y),gv=c(0,0,0.5,0.5,1,1,1,0.5,0.5,0))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS





