\name{cov.wml}
\alias{cov.wml}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Weighted Covariance Matrices (Maximum Likelihood)}
\description{
     Returns a list containing estimates of the weighted covariance
     matrix and the mean of the data, and optionally of the (weighted)
     correlation matrix. The
     covariance matrix is divided by the sum of the weights,
     corresponding to \code{n} and the ML-estimator in the case of equal
     weights, as opposed to \code{n-1} for \code{\link{cov.wt}}. 
}
\usage{
cov.wml(x, wt = rep(1/nrow(x), nrow(x)), cor = FALSE, center = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a matrix or data frame.  As usual, rows are observations and
          columns are variables.}
  \item{wt}{a non-negative and non-zero vector of weights for each
          observation.  Its length must equal the number of rows of
          \code{x}.}
  \item{cor}{A logical indicating whether the estimated correlation
          weighted matrix will be returned as well.}
  \item{center}{Either a logical or a numeric vector specifying the centers
          to be used when computing covariances.  If \code{TRUE}, the
          (weighted) mean of each variable is used, if `\code{FALSE}, zero is
          used.  If \code{center} is numeric, its length must equal the
          number of columns of \code{x}.}
}
\value{
  A list containing the following named components:
  \item{cov}{the estimated (weighted) covariance matrix.}
  \item{center}{an estimate for the center (mean) of the data.}
  \item{n.obs}{the number of observations (rows) in \code{x}.}
  \item{wt}{the weights used in the estimation.  Only returned if given
    as an argument.}
  \item{cor}{the estimated correlation matrix.  Only returned if `cor' is
          `TRUE'.}
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{\code{\link{cov.wt}}, \code{\link{cov}}, \code{\link{var}}}

\examples{
  x <- c(1,2,3,4,5,6,7,8,9,10)
  y <- c(1,2,3,8,7,6,5,8,9,10)
  cov.wml(cbind(x,y),wt=c(0,0,0,1,1,1,1,1,0,0))
  cov.wt(cbind(x,y),wt=c(0,0,0,1,1,1,1,1,0,0))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
