\name{piridge.zeroes}
\alias{piridge.zeroes}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extrema of two-component Gaussian mixture}
\description{
  By use of the Pi-function in Ray and Lindsay, 2005, locations of 
  two-component Gaussian mixture density extrema or saddlepoints are computed. 
}
\usage{
piridge.zeroes(prop, mu1, mu2, Sigma1, Sigma2, alphamin=0,
                          alphamax=1,by=0.001)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{prop}{proportion of mixture component 1.}
  \item{mu1}{mean vector of component 1.}
  \item{mu2}{mean vector of component 2.}
  \item{Sigma1}{covariance matrix of component 1.}
  \item{Sigma2}{covariance matrix of component 2.}
  \item{alphamin}{minimum alpha value.}
  \item{alphamax}{maximum alpha value.}
  \item{by}{interval between alpha-values where to look for extrema.}
}

\value{
  list with components
  \item{number.zeroes}{number of zeroes of Pi-function, i.e.,
    extrema or saddlepoints of density.}
  \item{estimated.roots}{estimated \code{alpha}-values at which extrema
  or saddlepoints occur.}
}

\references{
  Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  q <- piridge.zeroes(0.2,c(1,1),c(2,5),diag(2),diag(2),by=0.1)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



