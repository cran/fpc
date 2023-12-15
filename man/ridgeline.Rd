\name{ridgeline}
\alias{ridgeline}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ridgeline computation}
\description{
  Computes \eqn{(\alpha*\Sigma_1^{-1}+(1-\alpha)*\Sigma_2^{-1})^{-1}* 
    \alpha*(\Sigma_1^{-1}*\mu_1)+(1-\alpha)*(\Sigma_2^{-1}*\mu_2)}{%
  (alpha*Sigma1^{-1}+(1-alpha)*Sigma2^{-1})^{-1}* 
    alpha*(Sigma_1^{-1}*mu_1)+(1-alpha)*(Sigma_2^{-1}*mu_2)}
  as required for the
  computation of the ridgeline (Ray and Lindsay, 2005) to find
  all density extrema of a two-component Gaussian mixture with
  mean vectors mu1 and mu2 and covariance matrices Sigma1, Sigma2.
}
\usage{
ridgeline(alpha, mu1, mu2, Sigma1, Sigma2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{alpha}{numeric between 0 and 1.}
  \item{mu1}{mean vector of component 1.}
  \item{mu2}{mean vector of component 2.}
  \item{Sigma1}{covariance matrix of component 1.}
  \item{Sigma2}{covariance matrix of component 2.}
}

\value{
  A vector. See above.
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
  ridgeline(0.5,c(1,1),c(2,5),diag(2),diag(2))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



