\name{cmahal}
\alias{cmahal}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generation of tuning constant for Mahalanobis fixed point clusters.}
\description{
  Generates tuning constants \code{ca}
  for \code{\link{fixmahal}} dependent on
  the number of points and variables of the current fixed point cluster
  (FPC).
  
  This is experimental and only thought for use in \code{\link{fixmahal}}.
}
\usage{
cmahal(n, p, nmin, cmin, nc1, c1 = cmin, q = 1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Number of points.}
  \item{p}{positive integer. Number of variables.}
  \item{nmin}{integer larger than 1. Smallest number of points for which
    \code{ca} is computed. For smaller FPC sizes, \code{ca} is set to
    the value for \code{nmin}.}
  \item{cmin}{positive number. Minimum value for \code{ca}.}
  \item{nc1}{positive integer. Number of points at which \code{ca=c1}.}
  \item{c1}{positive numeric. Tuning constant for \code{cmahal}.
    Value for \code{ca} for FPC size equal to \code{nc1}.}
  \item{q}{numeric between 0 and 1. 1 for steepest possible descent of
    \code{ca} as function of the FPC size. Should presumably always be 1.}
}
\details{
  Some experiments suggest that the tuning constant \code{ca} should
  decrease with increasing FPC size and increase with increasing
  \code{p} in \code{\link{fixmahal}}. This is to prevent too small
  meaningless FPCs while maintaining the significant larger
  ones. \code{cmahal} with \code{q=1} computes \code{ca} in such a way
  that as long as \code{ca>cmin}, the decrease in \code{n} is as steep
  as possible in order to maintain the validity of the convergence
  theorem in Hennig and Christlieb (2002).
}
\value{
  A numeric vector of length \code{n}, giving the values for \code{ca}
  for all FPC sizes smaller or equal to \code{n}.
}
\references{
  Hennig, C. and Christlieb, N. (2002) Validating visual clusters in
  large datasets: Fixed point clusters of spectral features,
  \emph{Computational Statistics and Data Analysis} 40, 723-739.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{\code{\link{fixmahal}}}

\examples{
  plot(1:100,cmahal(100,3,nmin=5,cmin=qchisq(0.99,3),nc1=90),
       xlab="FPC size", ylab="cmahal")
}
\keyword{cluster}% at least one, from doc/KEYWORDS

