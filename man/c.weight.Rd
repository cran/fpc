\name{c.weight}
\alias{c.weight}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Weight function for AWC}
\description{
  For use in \code{awcoord} only.
}
\usage{
c.weight(x, ca) 

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
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{awcoord}}
}

\examples{
  c.weight(4,1)
}
\keyword{arith}% at least one, from doc/KEYWORDS



