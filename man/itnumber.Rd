\name{itnumber}
\alias{itnumber}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Number of regression fixed point cluster iterations}
\description{
  Computes the number of fixed point iterations needed by
  \code{\link{fixreg}} to find \code{mtf} times
  a fixed point cluster (FPC) of size
  \code{cn} with an approximated probability of \code{prob}.

  Thought for use within \code{\link{fixreg}}.
}
\usage{
itnumber(n, p, cn, mtf, prob = 0.95, maxir = 20000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Total number of points.}
  \item{p}{positive integer. Number of independent variables.}
  \item{cn}{positive integer smaller or equal to \code{n}.
    Size of the FPC.}
  \item{mtf}{positive integer.}
  \item{prob}{number between 0 and 1.}
  \item{maxir}{positive integer. \code{itnumber} is set to this value if
    it would otherwise be larger.}
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
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{\code{\link{fixreg}}, \code{\link{clusexpect}}}

\examples{
  itnumber(500,4,150,2)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
