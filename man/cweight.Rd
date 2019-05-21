\name{cweight}
\alias{cweight}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Weight function for AWC}
\description{
  For use in \code{awcoord} only.
}
\usage{
cweight(x, ca) 

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{numerical.}
  \item{ca}{numerical.}
}
% \details{
% }
\value{
  \code{ca/x} if smaller than 1, else 1.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{awcoord}}
}

\examples{
  cweight(4,1)
}
\keyword{arith}% at least one, from doc/KEYWORDS



