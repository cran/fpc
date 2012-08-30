\name{clusexpect}
\alias{clusexpect}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Expected value of the number of times a fixed point
  cluster is found}
\description{
  A rough approximation of the expectation of the number of times a well
  separated fixed point
  cluster (FPC) of size \code{n} is found in \code{ir} fixed point
  iterations of \code{\link{fixreg}}.
}
\usage{
clusexpect(n, p, cn, ir)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Total number of points.}
  \item{p}{positive integer. Number of independent variables.}
  \item{cn}{positive integer smaller or equal to \code{n}.
    Size of the FPC.}
  \item{ir}{positive integer. Number of fixed point iterations.}
}
\details{
  The approximation is based on the assumption that a well separated FPC
  is found iff all \code{p+2} points of the initial coinfiguration come
  from the FPC. The value is \code{ir} times the probability for
  this. For a discussion of this assumption cf. Hennig (2002).
}
\value{
  A number.
}

\references{
  Hennig, C. (2002) Fixed point clusters for linear regression:
  computation and comparison, \emph{Journal of
    Classification} 19, 249-276.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{\code{\link{fixreg}}}

\examples{
  clusexpect(500,4,150,2000)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{cluster}
