\name{bhattacharyya.dist}
\alias{bhattacharyya.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Bhattacharyya distance between Gaussian distributions}
\description{
  Computes Bhattacharyya distance between two multivariate
  Gaussian distributions. See Fukunaga (1990).
}
\usage{
bhattacharyya.dist(mu1, mu2, Sigma1, Sigma2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mu1}{mean vector of component 1.}
  \item{mu2}{mean vector of component 2.}
  \item{Sigma1}{covariance matrix of component 1.}
  \item{Sigma2}{covariance matrix of component 2.}
}

\value{
  The Bhattacharyya distance between the two Gaussian distributions.
}

\references{
Fukunaga, K. (1990) \emph{Introduction to Statistical Pattern
Recognition}, 2nd edition, Academic
Press, New York.

  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  round(bhattacharyya.dist(c(1,1),c(2,5),diag(2),diag(2)),digits=2)
}
\keyword{multivariate}



