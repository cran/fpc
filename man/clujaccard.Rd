\name{clujaccard}
\alias{clujaccard}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Jaccard similarity between logical vectors}
\description{
  Jaccard similarity between logical or 0-1 vectors:
  \code{sum(c1 & c2)/sum(c1 | c2)}.
}
\usage{
clujaccard(c1,c2,zerobyzero=NA)
}
\arguments{
  \item{c1}{logical or 0-1-vector.}
  \item{c2}{logical or 0-1-vector (same length).}
  \item{zerobyzero}{result if \code{sum(c1 | c2)=0}.}
}

\value{
  Numeric between 0 and 1.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  c1 <- rep(TRUE,10)
  c2 <- c(FALSE,rep(TRUE,9))
  clujaccard(c1,c2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS



