\name{unimodal.ind}
\alias{unimodal.ind}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Is a fitted denisity unimodal or not?}
\description{
  Checks whether a series of fitted density values (such as given out as
  \code{y}-component of  \code{\link{density}}) is unimodal.
}
\usage{
  unimodal.ind(y)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{numeric vector of fitted density values in order of
    increasing x-values such as given out as
  \code{y}-component of  \code{\link{density}}.}
}

\value{
  Logical. \code{TRUE} if unimodal.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
unimodal.ind(c(1,3,3,4,2,1,0,0))
}
\keyword{univar}% at least one, from doc/KEYWORDS
% \keyword{multivariate}



