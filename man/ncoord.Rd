\name{ncoord}
\alias{ncoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Neighborhood based discriminant coordinates}
\description{
  Neighborhood based discriminant coordinates as defined in Hastie and
  Tibshirani (1996) and a robustified version as defined in Hennig (2003).
  The principle is to maximize the projection of a between
  classes covariance matrix, which is defined by averaging the
  between classes covariance matrices in the neighborhoods of all points.  
}
\usage{
ncoord(xd, clvecd, nn=50, weighted=FALSE,
                    sphere="mcd", orderall=TRUE, countmode=1000, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{nn}{integer. Number of points which belong to the neighborhood
    of each point (including the point itself).}
  \item{weighted}{logical. \code{FALSE} corresponds to the original
    method of Hastie and Tibshirani (1996). If \code{TRUE},
    the between classes
    covariance matrices B are weighted by w/trace B, where w is some
    weight depending on the sizes of the
    classes in the neighborhood. Division by trace B reduces the effect
    of outliers. \code{TRUE} cooresponds to WNC as defined in Hennig
    (2003).}
  \item{sphere}{a covariance matrix or one of
    "mve", "mcd", "classical", "none". The matrix used for sphering the
    data. "mcd" and "mve" are robust covariance matrices as implemented
    in \code{\link{cov.rob}}. "classical" refers to the classical
    covariance matrix. "none" means no sphering and use of the raw
    data.}
  \item{orderall}{logical. By default, the neighborhoods are computed by
    ordering all points each time. If \code{FALSE}, the neighborhoods
    are computed by selecting \code{nn} times the nearest point from the
    remaining points, which may be faster sometimes.}
  \item{countmode}{optional positive integer. Every \code{countmode}
    algorithm runs \code{ncoord} shows a message.}
  \item{...}{no effect}
}
% \details{
% }
\value{
  List with the following components
  \item{ev}{eigenvalues in descending order.}
  \item{units}{columns are coordinates of projection basis vectors.
  New points \code{x} can be projected onto the projection basis vectors
  by \code{x \%*\% units}}
  \item{proj}{projections of \code{xd} onto \code{units}.} 
}
\references{
  Hastie, T. and Tibshirani, R.  (1996). Discriminant adaptive nearest
  neighbor classification. \emph{IEEE Transactions on Pattern Analysis
  and Machine Intelligence} 18, 607-616. 

  Hennig, C. (2003) Symmetric, asymmetric, and robust linear dimension
  reduction for classification, submitted,
  \url{http://stat.ethz.ch/Research-Reports/108.html}.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}
}

\seealso{
  \code{\link{plotcluster}} for straight forward discriminant plots.
  \code{\link{discrproj}} for alternatives.
  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
  set.seed(4634)
  face <- rFace(600,dMoNo=2,dNoEy=0)
  grface <- as.integer(attr(face,"grouping"))
  ncf <- ncoord(face,grface)
  plot(ncf$proj,col=grface)
  ncf2 <- ncoord(face,grface,weighted=TRUE)
  plot(ncf2$proj,col=grface)
  # ...done in one step by function plotcluster.
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



