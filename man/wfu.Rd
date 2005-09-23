\name{wfu}
\alias{wfu}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Weight function (for Mahalabobis distances)}
\description{
  Function of the elements of \code{md}, which is 1 for arguments smaller
  than \code{ca}, 0 for arguments larger than \code{ca2} and linear
  (default: continuous) in between.
  
  Thought for use in \code{fixmahal}.
}
\usage{
wfu(md, ca, ca2, a1 = 1/(ca - ca2), a0 = -a1 * ca2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{md}{vector of positive numericals.}
  \item{ca}{positive numerical.}
  \item{ca2}{positive numerical.}
  \item{a1}{numerical. Slope.}
  \item{a0}{numerical. Intercept.}
}

\value{
  A vector of numericals between 0 and 1.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{\code{\link{fixmahal}}}

\examples{
  md <- seq(0,10,by=0.1)
  wfu(md,ca=5,ca2=8)
}
\keyword{arith}% at least one, from doc/KEYWORDS







