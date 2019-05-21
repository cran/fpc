\name{mahalconf}
\alias{mahalconf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mahalanobis fixed point clusters initial configuration}
\description{
  Generates an initial configuration of \code{startn} points from
  dataset \code{x} for the \code{\link{fixmahal}}
  fixed point iteration.

  Thought only for use within \code{\link{fixmahal}}.
}
\usage{
mahalconf(x, no, startn, covall, plot)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{numerical matrix. Rows are points, columns are variables.}
  \item{no}{integer between 1 and \code{nrow(x)}. Number of the first
    point of the configuration.}
  \item{startn}{integer between 1 and \code{nrow(x)}.}
  \item{covall}{covariance matrix for the computation of the first
    Mahalanobis distances.}
  \item{plot}{a string. If equal to \code{"start"} or \code{"both"},the
    first two variables and the first \code{ncol(x)+1} points are plotted.}
}
\details{
  \code{mahalconf} first chooses the \eqn{p} (number of variables)
  nearest points to point no. \code{no} in terms of the Mahalanobis
  distance w.r.t. \code{covall}, so that there are \eqn{p+1} points.
  In every further step, the covariance
  matrix of the current configuration is computed and the nearest point
  in terms of the new Mahalanobis distance is
  added. \code{\link{solvecov}} is used to invert singular covariance
  matrices. 
}
\value{
  A logical vector of length \code{nrow(x)}.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}}

\seealso{\code{\link{fixmahal}}, \code{\link{solvecov}}}

\examples{
  set.seed(4634)
  face <- rFace(600,dMoNo=2,dNoEy=0,p=2)
  mahalconf(face,no=200,startn=20,covall=cov(face),plot="start")
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}
