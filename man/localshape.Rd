\name{localshape}
\alias{localshape}
\title{Local shape matrix}
\description{
  This computes a matrix formalising 'local shape', i.e., aggregated
  standardised variance/covariance in a Mahalanobis neighbourhood of the data
  points. This can be used for finding clusters when used as one of the
  covariance matrices in 
  Invariant Coordinate Selection (function \code{ics} in package
  \code{ICS}), see Hennig's
  discussion and rejoinder of Tyler et al. (2009).
}
\usage{
  localshape(xdata,proportion=0.1,mscatter="mcd",mcdalpha=0.8,
                       covstandard="det")
}
\arguments{
  \item{xdata}{objects times variables data matrix.}
  \item{proportion}{proportion of points to be considered as neighbourhood.}
  \item{mscatter}{"mcd" or "cov"; specified minimum covariance
    determinant or
    classical covariance matrix to be used for Mahalanobis distance
    computation.}
  \item{mcdalpha}{if \code{mscatter="mcd"}, this is the alpha parameter
    to be used by the MCD covariance matrix, i.e. one minus the
    asymptotic breakdown point, see \code{\link[robustbase]{covMcd}}.}
  \item{covstandard}{one of "trace", "det" or "none", determining by
    what constant the pointwise neighbourhood covariance matrices are
    standardised. "det" makes the affine equivariant, as noted in the
    discussion rejoinder of Tyler et al. (2009).}
}

\value{
  The local shape matrix.
}

\references{
  Tyler, D. E., Critchley, F., Duembgen, L., Oja, H. (2009)
  Invariant coordinate selection (with discussion).
  \emph{Journal of the Royal Statistical Society, Series B}, 549-592.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}
\examples{
  options(digits=3)
  data(iris)
  localshape(iris[,-5],mscatter="cov")
}
\keyword{multivariate}

