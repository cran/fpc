\name{piridge}
\alias{piridge}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ridgeline Pi-function}
\description{
  The Pi-function is given in (6) in Ray and Lindsay, 2005. Equating it
  to the mixture proportion yields locations of two-component Gaussian
  mixture density extrema. 
}
\usage{
piridge(alpha, mu1, mu2, Sigma1, Sigma2, showplot=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{alpha}{sequence of values between 0 and 1 for which the Pi-function
is computed.}
  \item{mu1}{mean vector of component 1.}
  \item{mu2}{mean vector of component 2.}
  \item{Sigma1}{covariance matrix of component 1.}
  \item{Sigma2}{covariance matrix of component 2.}
  \item{showplot}{logical. If \code{TRUE}, the Pi-function is plotted
against \code{alpha}.}
}

\value{
  Vector of values of the Pi-function for values of \code{alpha}.
}

\references{
  Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
  q <- piridge(seq(0,1,0.1),c(1,1),c(2,5),diag(2),diag(2))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



