\name{pamk}
\alias{pamk}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Partitioning around medoids with estimation of number of clusters}
\description{
  This calls the function \code{\link[cluster]{pam}} or
  \code{\link[cluster]{clara}} to perform a
  partitioning around medoids clustering with the number of clusters
  estimated by optimum average silhouette width (see
  \code{\link[cluster]{pam.object}}) or Calinski-Harabasz
  index (\code{\link{calinhara}}). The Duda-Hart test
  (\code{\link{dudahart2}}) is applied to decide whether there should be
  more than one cluster (unless 1 is excluded as number of clusters or
  data are dissimilarities).  
}
\usage{
pamk(data,krange=2:10,criterion="asw", usepam=TRUE,
     scaling=FALSE, alpha=0.001, diss=inherits(data, "dist"),
     critout=FALSE, ns=10, seed=NULL, ...)
}
\arguments{
  \item{data}{a data matrix or data frame or something that can be
    coerced into a matrix, or dissimilarity matrix or
    object. See \code{\link[cluster]{pam}} for more information.}
  \item{krange}{integer vector. Numbers of clusters which are to be
    compared by the average silhouette width criterion. Note: average
    silhouette width and Calinski-Harabasz can't estimate number of
    clusters \code{nc=1}. If 1 is included, a Duda-Hart test is applied
    and 1 is estimated if this is not significant.}
  \item{criterion}{one of \code{"asw"}, \code{"multiasw"} or
    \code{"ch"}. Determines whether average silhouette width (as given
    out by  \code{\link[cluster]{pam}}/\code{\link[cluster]{clara}}, or
    as computed by \code{\link{distcritmulti}} if \code{"multiasw"} is
    specified; recommended for large data sets with \code{usepam=FALSE})
    or Calinski-Harabasz is applied. Note that the original
    Calinski-Harabasz index is not defined for dissimilarities; if
    dissimilarity data is run with \code{criterion="ch"}, the
    dissimilarity-based generalisation in Hennig and Liao (2013) is
    used.}
  \item{usepam}{logical. If \code{TRUE}, \code{\link[cluster]{pam}} is
    used, otherwise \code{\link[cluster]{clara}} (recommended for large
    datasets with 2,000 or more observations; dissimilarity matrices can
    not be used with \code{\link[cluster]{clara}}).}
  \item{scaling}{either a logical value or a numeric vector of length
    equal to the number of variables. If \code{scaling} is a numeric
    vector with length equal to the number of variables, then each
    variable is divided by the corresponding value from \code{scaling}.
    If \code{scaling} is \code{TRUE} then scaling is done by dividing
    the (centered) variables by their root-mean-square, and if
    \code{scaling} is \code{FALSE}, no scaling is done.}
  \item{alpha}{numeric between 0 and 1, tuning constant for
    \code{\link{dudahart2}} (only used for 1-cluster test).} 
  \item{diss}{logical flag: if \code{TRUE} (default for \code{dist} or
    \code{dissimilarity}-objects), then \code{data} will be considered
    as a dissimilarity matrix (and the potential number of clusters 1
    will be ignored).  If \code{FALSE}, then \code{data} will
    be considered as a matrix of observations by variables.}
  \item{critout}{logical. If \code{TRUE}, the criterion value is printed
    out for every number of clusters.}
  \item{ns}{passed on to \code{\link{distcritmulti}} if
    \code{criterion="multiasw"}.}
  \item{seed}{passed on to \code{\link{distcritmulti}} if
    \code{criterion="multiasw"}.}
  \item{...}{further arguments to be transferred to
    \code{\link[cluster]{pam}} or \code{\link[cluster]{clara}}.}
}

\note{
  \code{\link[cluster]{clara}} and \code{\link[cluster]{pam}}
  can handle \code{NA}-entries (see their documentation) but
  \code{\link{dudahart2}} cannot. Therefore \code{NA} should not occur
  if 1 is in \code{krange}.
}
  
\value{
  A list with components
  \item{pamobject}{The output of the optimal run of the
    \code{\link[cluster]{pam}}-function.}
  \item{nc}{the optimal number of clusters.}
  \item{crit}{vector of criterion values for numbers of
    clusters. \code{crit[1]} is the p-value of the Duda-Hart test
    if 1 is in \code{krange} and \code{diss=FALSE}.}
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\references{
  Calinski, R. B., and Harabasz, J. (1974) A Dendrite Method for Cluster 
Analysis, \emph{Communications in Statistics}, 3, 1-27.

  Duda, R. O. and Hart, P. E. (1973) \emph{Pattern Classification and
  Scene Analysis}. Wiley, New York.

  Hennig, C. and Liao, T. (2013) How to find an appropriate clustering
  for mixed-type variables with application to socio-economic
  stratification, \emph{Journal of the Royal Statistical Society, Series
  C Applied Statistics}, 62, 309-369.

  Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
     An Introduction to Cluster Analysis". Wiley, New York.
}
\seealso{
  \code{\link[cluster]{pam}}, \code{\link[cluster]{clara}}
  \code{\link{distcritmulti}}
}
\examples{
  options(digits=3)
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  pk1 <- pamk(face,krange=1:5,criterion="asw",critout=TRUE)
  pk2 <- pamk(face,krange=1:5,criterion="multiasw",ns=2,critout=TRUE)
# "multiasw" is better for larger data sets, use larger ns then.
  pk3 <- pamk(face,krange=1:5,criterion="ch",critout=TRUE)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



