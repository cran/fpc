\name{jittervar}
\alias{jittervar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Jitter variables in a data matrix}
\description{
  Jitters some variables in a data matrix.
}
\usage{
  jittervar(x,jitterv=NULL,factor=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame.}
  \item{jitterv}{vector of numbers of variables to be jittered.}
  \item{factor}{numeric. Passed on to \code{\link{jitter}}. See the
    documentation there. The higher, the more jittering.}
}

\value{
  data matrix or data frame with jittered variables.
  }

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{jitter}}}

\examples{
  set.seed(776655)
  v1 <- rnorm(20)
  v2 <- rnorm(20)
  d1 <- sample(1:5,20,replace=TRUE)
  d2 <- sample(1:4,20,replace=TRUE)
  ldata <- cbind(v1,v2,d1,d2)
  jv <- jittervar(ldata,jitterv=3:4)
}

\keyword{manip}% __ONLY ONE__ keyword per line
