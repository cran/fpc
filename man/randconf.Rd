\name{randconf}
\alias{randconf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate a sample indicator vector}
\description{
  Generates a logical vector of length \code{n} with \code{p TRUE}s.
}
\usage{
randconf(n, p)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer.}
  \item{p}{positive integer.}
}

\value{
  A logical vector.
}

\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{\code{\link{sample}}}

\examples{
  randconf(10,3)
}
\keyword{distribution}% at least one, from doc/KEYWORDS

