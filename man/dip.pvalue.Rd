\name{dip.pvalue}
\alias{dip.pvalue}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Approximates p-value for dip test}
\description{
  Approximates p-value for dip test (see \code{\link[diptest]{dip}})
  using the information in \code{\link[diptest]{qDiptab}}.
}
\usage{
  dip.pvalue(d,n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{numeric. Value of dip statistic.}
  \item{n}{integer. Number of observations.}
}

\value{
  Approximated p-value.
}

\references{
     J. A. Hartigan and P. M. Hartigan (1985) The Dip Test of
     Unimodality, \emph{Annals of Statistics}, 13, 70-84.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  dip.pvalue(0.08,15)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
% \keyword{multivariate}



