\name{tdecomp}
\alias{tdecomp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Root of singularity-corrected eigenvalue decomposition}
\description{
  Computes transposed eigenvectors of matrix \code{m} times diagonal of
  square root of eigenvalues so that eigenvalues smaller than 1e-6 are
  set to 1e-6.
}
\usage{
tdecomp(m)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{m}{a symmetric matrix of minimum format 2*2.}
}
\details{
  Thought for use in \code{discrcoord} only.}
\value{
  a matrix.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\note{
  Thought for use within \code{\link{discrcoord}} only.
}

\examples{
x <- rnorm(10)
y <- rnorm(10)
z <- cov(cbind(x,y))
tdecomp(z)
}
\keyword{array}% at least one, from doc/KEYWORDS

