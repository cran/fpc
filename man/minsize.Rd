\name{minsize}
\alias{minsize}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Minimum size of regression fixed point cluster}
\description{
  Computes the minimum size of a fixed point cluster (FPC) which is
  found at least \code{mtf} times with approximated
  probability \code{prob} by
  \code{ir} fixed point iterations of \code{\link{fixreg}}.

  Thought for use within \code{\link{fixreg}}.
}
\usage{
minsize(n, p, ir, mtf, prob = 0.5)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Total number of points.}
  \item{p}{positive integer. Number of independent variables.}
  \item{ir}{positive integer. Number of fixed point iterations.}
  \item{mtf}{positive integer.}
  \item{prob}{numerical between 0 and 1.}
}
\details{
  The computation is based on the binomial distribution with probability
  given by \code{\link{clusexpect}} with \code{ir=1}.
}
\value{
  An integer.
}
\references{
  Hennig, C. (2002) Fixed point clusters for linear regression:
  computation and comparison, \emph{Journal of
    Classification} 19, 249-276.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{\code{\link{fixreg}}, \code{\link{clusexpect}},
  \code{\link{itnumber}}}

\examples{
  minsize(500,4,7000,2)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
