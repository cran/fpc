\name{mahalanodisc}
\alias{mahalanodisc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mahalanobis for AWC}
\description{
  Vector of Mahalanobis distances or their root. For use in \code{awcoord} only.
}
\usage{
mahalanodisc(x2, mg, covg, modus="square") 

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x2}{numerical data matrix.}
  \item{mg}{mean vector.}
  \item{covg}{covariance matrix.}
  \item{modus}{"md" (roots of Mahalanobis distances) or "square"
    (original squared form of Mahalanobis distances).}
}
\details{
  The covariance matrix
  is inverted by use of
  \code{\link{solvecov}}, which can be expected to give
  reasonable results for singular within-class covariance matrices.
}
% \details{
% }
\value{
  vector of (rooted) Mahalanobis distances.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{awcoord}}, \code{\link{solvecov}}
}

\examples{
  x <- cbind(rnorm(50),rnorm(50))
  mahalanodisc(x,c(0,0),cov(x))
  mahalanodisc(x,c(0,0),matrix(0,ncol=2,nrow=2))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS



