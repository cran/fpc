\name{bhattacharyya.matrix}
\alias{bhattacharyya.matrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Matrix of pairwise Bhattacharyya distances}
\description{
  Computes Bhattachryya distances for pairs of components
  given the parameters of a Gaussian mixture. 
}
\usage{
bhattacharyya.matrix(muarray,Sigmaarray,ipairs="all", 
                                 misclassification.bound=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{muarray}{matrix of component means (different components are in
    different columns).}
  \item{Sigmaarray}{three dimensional array with component covariance
    matrices (the third dimension refers to components).}
  \item{ipairs}{\code{"all"} or list of vectors of two integers. If
    \code{ipairs="all"}, computations are carried out for all pairs of
    components. Otherwise, ipairs gives the pairs of components for
    which computations are carried out.}
  \item{misclassification.bound}{logical. If \code{TRUE}, upper bounds
    for misclassification probabilities \code{exp(-b)}
    are given out instead of the original Bhattacharyya distances \code{b}.}  
}

\value{
  A matrix with Bhattacharyya distances (or derived misclassification
  bounds, see above) between pairs of Gaussian distributions with the
  provided parameters. If \code{ipairs!="all"}, the Bhattacharyya
  distance is taken as 0 and the misclassification bound as 1.
}

\references{
  Fukunaga, K. (1990) \emph{Introduction to Statistical Pattern
Recognition}, 2nd edition, Academic
Press, New York.
 
  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
  }
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{bhattacharyya.dist}}
  }

\examples{
  muarray <- cbind(c(0,0),c(0,0.1),c(10,10))
  sigmaarray <- array(c(diag(2),diag(2),diag(2)),dim=c(2,2,3))
  bhattacharyya.matrix(muarray,sigmaarray,ipairs=list(c(1,2),c(2,3)))

}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



