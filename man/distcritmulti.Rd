\name{distcritmulti}
\alias{distcritmulti}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distance based validity criteria for large data sets}
\description{
  Approximates average silhouette width or the Pearson version of
  Hubert's gamma criterion by hacking the
  dataset into pieces and averaging the subset-wise values, see Hennig
  and Liao (2013).
}  
\usage{
distcritmulti(x,clustering,part=NULL,ns=10,criterion="asw",
                    fun="dist",metric="euclidean",
                     count=FALSE,seed=NULL,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{cases times variables data matrix.}
  \item{clustering}{vector of integers indicating the clustering.}
  \item{part}{vector of integer subset sizes; sum should be smaller or
     equal to the number of cases of \code{x}. If \code{NULL}, subset sizes are
     chosen approximately equal.} 
  \item{ns}{integer. Number of subsets, only used if \code{part==NULL}.}
  \item{criterion}{\code{"asw"} or \code{"pearsongamma"}, specifies
      whether the average silhouette width or the Pearson version of
      Hubert's gamma is computed.}
  \item{fun}{\code{"dist"} or \code{"daisy"}, specifies
      which function is used for computing dissimilarities.}
  \item{metric}{passed on to \code{\link{dist}} (as argument
      \code{method}) or \code{\link[cluster]{daisy}} to determine which
    dissimilarity is used.}
  \item{count}{logical. if \code{TRUE}, the subset number just processed
    is printed.}
  \item{seed}{integer, random seed. (If \code{NULL}, result depends on
    random numbers.)}
  \item{...}{further arguments to be passed on to \code{\link{dist}} or
    \code{\link[cluster]{daisy}}.}
}

\value{
  A list with components \code{crit.overall,crit.sub,crit.sd,part}.
  \item{crit.overall}{value of criterion.}
  \item{crit.sub}{vector of subset-wise criterion values.}
  \item{crit.sd}{standard deviation of \code{crit.sub}, can be used to
    assess stability.}
  \item{subsets}{list of case indexes in subsets.}
}


\references{
  Halkidi, M., Batistakis, Y., Vazirgiannis, M. (2001) On Clustering
  Validation Techniques, \emph{Journal of Intelligent Information
    Systems}, 17, 107-145.
  
  Hennig, C. and Liao, T. (2013) How to find an appropriate clustering
  for mixed-type variables with application to socio-economic
  stratification, \emph{Journal of the Royal Statistical Society, Series
  C Applied Statistics}, 62, 309-369.

    Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
  An Introduction to Cluster Analysis". Wiley, New York.
  }

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{cluster.stats}}, \code{\link[cluster]{silhouette}}}

\examples{
    set.seed(20000)
    options(digits=3)
    face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
    clustering <- as.integer(attr(face,"grouping"))
    distcritmulti(face,clustering,ns=3,seed=100000,criterion="pearsongamma")
}

\keyword{cluster}% __ONLY ONE__ keyword per line
