\name{zmisclassification.matrix}
\alias{zmisclassification.matrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Matrix of misclassification probabilities between mixture components}
\description{
   Matrix of misclassification probabilities in a mixture distribution
   between two mixture components from estimated posterior probabilities
   regardless of component parameters, see Hennig (2010).
}
\usage{
zmisclassification.matrix(z,pro=NULL,clustering=NULL,
                                      ipairs="all",symmetric=TRUE,
                                      stat="max")

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{z}{matrix of posterior probabilities for observations (rows) to
    belong to mixture components (columns), so entries need to sum up to
    1 for each row.}
  \item{pro}{vector of component proportions, need to sum up to
    1. Computed from \code{z} as default.}
  \item{clustering}{vector of integers giving the estimated mixture
    components for every observation. Computed from \code{z} as
    default.}
  \item{ipairs}{\code{"all"} or list of vectors of two integers. If
    \code{ipairs="all"}, computations are carried out for all pairs of
    components. Otherwise, ipairs gives the pairs of components for
    which computations are carried out.}
  \item{symmetric}{logical. If \code{TRUE}, the matrix is symmetrised,
    see parameter \code{stat}.}
  \item{stat}{\code{"max"} or \code{"mean"}. The statistic by which the
    two misclassification probabilities are aggregated if
    \code{symmetric=TRUE}.} 
}

\value{
  A matrix with the (symmetrised, if required) misclassification
  probabilities between each pair of mixture components. If
  \code{symmetric=FALSE}, matrix entry \code{[i,j]} is the estimated
  probability that an observation generated by component
  \code{j} is classified to component \code{i} by maximum a posteriori rule.
}

\references{
  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{confusion}}
  }
\examples{
  set.seed(12345)
  m <- rpois(20,lambda=5)
  dim(m) <- c(5,4)
  m <- m/apply(m,1,sum)
  round(zmisclassification.matrix(m,symmetric=FALSE),digits=2) 
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



