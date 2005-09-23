\name{can}
\alias{can}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generation of the tuning constant for regression fixed point clusters}
\description{
  Generates tuning constants \code{ca}
  for \code{\link{fixreg}} dependent on
  the number of points and variables of the dataset.
  
  Only thought for use in \code{\link{fixreg}}.
  
}
\usage{
can(n, p)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Number of points.}
  \item{p}{positive integer. Number of independent variables.}
}
\details{
  The formula is
  \eqn{3+33/(n*2^{-(p-1)/2})^{1/3}+2900000/(n*2^{-(p-1)/2})^3}. For
  justification cf. Hennig (2002).
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
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{\code{\link{fixreg}}}

\examples{
  can(429,3)
}
\keyword{arith}% at least one, from doc/KEYWORDS

