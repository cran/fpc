\name{dridgeline}
\alias{dridgeline}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density along the ridgeline}
\description{
  Computes the density of a two-component Gaussian mixture along the 
  ridgeline (Ray and Lindsay, 2005), along which 
  all its density extrema are located.
}
\usage{
dridgeline(alpha=seq(0,1,0.001), prop,
                          mu1, mu2, Sigma1, Sigma2, showplot=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{alpha}{sequence of values between 0 and 1 for which the density
is computed.}
  \item{prop}{mixture proportion of first component.}
  \item{mu1}{mean vector of component 1.}
  \item{mu2}{mean vector of component 2.}
  \item{Sigma1}{covariance matrix of component 1.}
  \item{Sigma2}{covariance matrix of component 2.}
  \item{showplot}{logical. If \code{TRUE}, the density is plotted
against \code{alpha}.}
  \item{...}{further arguments to be passed on to plot.}
}

\value{
  Vector of density values for values of \code{alpha}.
}

\references{
  Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  q <- dridgeline(seq(0,1,0.1),0.5,c(1,1),c(2,5),diag(2),diag(2))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



