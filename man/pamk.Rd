\name{pamk}
\alias{pamk}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Partitioning around medoids with estimation of number of clusters}
\description{
  This calls the function \code{\link[cluster]{pam}} to perform a
  partitioning around medoids clustering with the number of clusters
  estimated by optimum average silhouette width.
}
\usage{
pamk(data,krange=2:10,scaling=FALSE, diss=inherits(data, "dist"),...)
}
\arguments{
  \item{data}{a data matrix or data frame, or dissimilarity matrix or
    object. See \code{\link[cluster]{pam}} for more information.}
  \item{krange}{integer vector. Numbers of clusters which are to be
    compared by the average silhouette width criterion. Note: This can't
    estimate number of clusters \code{nc=1}, and therefore 1 should not be in
    \code{krange}.}
  \item{scaling}{either a logical value or a numeric vector of length
    equal to the number of variables. If \code{scaling} is a numeric
    vector with length equal to the number of variables, then each
    variable is divided by the corresponding value from \code{scaling}.
    If \code{scaling} is \code{TRUE} then scaling is done by dividing
    the (centered) variables by their root-mean-square, and if
    \code{scaling} is \code{FALSE}, no scaling is done.}
  \item{diss}{logical flag: if \code{TRUE} (default for \code{dist} or
    \code{dissimilarity}-objects), then \code{data} will be considered
    as a dissimilarity matrix.  If \code{FALSE}, then \code{data} will
    be considered as a matrix of observations by variables.}
  \item{...}{further arguments to be transferred to
    \code{\link[cluster]{pam}}.}
}

\value{
  A list with components
  \item{pamobject}{The output of the optimal run of the
    \code{\link[cluster]{pam}}-function.}
  \item{nc}{the optimal number of clusters.}
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\references{
     Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
     An Introduction to Cluster Analysis". Wiley, New York.
}
\seealso{
  \code{\link[cluster]{pam}}
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  pk <- pamk(face)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



