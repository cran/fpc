\name{sseg}
\alias{sseg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Position in a similarity vector}
\description{
  \code{sseg(i,j)} gives the position of the similarity of objects
  \code{i} and \code{j} in the similarity vectors produced by
  \code{fixreg} and \code{fixmahal}.
  \code{sseg} should only be used as an auxiliary function in
  \code{fixreg} and \code{fixmahal}.
}
\usage{
sseg(i, j)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{i}{positive integer.}
  \item{j}{positive integer.}
}
\value{A positive integer.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\examples{
sseg(3,4)
}
\keyword{utilities}% at least one, from doc/KEYWORDS
