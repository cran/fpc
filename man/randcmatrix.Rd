\name{randcmatrix}
\alias{randcmatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Random partition matrix}
\description{
  For use within \code{regmix}. Generates a random
  0-1-matrix with \code{n} rows
  and \code{cln} columns so that every row contains exactly one one and
  every columns contains at least \code{p+3} ones.
}
\usage{
randcmatrix(n,cln,p)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{positive integer. Number of rows.}
  \item{cln}{positive integer. Number of columns.}
  \item{p}{positive integer. See above.}
}
\value{
  An \code{n*cln}-matrix.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}
\seealso{
  \code{\link{regmix}}
}
\examples{
set.seed(111)
randcmatrix(10,2,1)
}
\keyword{cluster}% at least one, from doc/KEYWORDS

