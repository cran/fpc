\name{diptest.multi}
\alias{diptest.multi}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Diptest for discriminant coordinate projection}
\description{
  Diptest (Hartigan and Hartigan, 1985, see \code{\link[diptest]{dip}})
  for data projected in discriminant coordinate separating optimally two
  class means (see \code{discrcoord}) as suggested by Tantrum, Murua and
  Stuetzle (2003).
}
\usage{
  diptest.multi(xdata,class,pvalue="uniform",M=100)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{matrix. Potentially multidimensional dataset.}
  \item{class}{vector of integers giving class numbers for observations.}
  \item{pvalue}{\code{"uniform"} or \code{"tantrum"}. Defines whether
    the p-value is computed from a uniform null model as suggested in
    Hartigan and Hartigan (1985, using \code{\link[diptest]{dip.test}}) or as
    suggested in Tantrum et al. (2003, using \code{dipp.tantrum}).}
  \item{M}{integer. Number of artificial datasets generated in order to
    estimate the p-value if \code{pvalue="tantrum"}.}
}

\value{
  The resulting p-value.
}

\references{
     J. A. Hartigan and P. M. Hartigan (1985) The Dip Test of
     Unimodality, \emph{Annals of Statistics}, 13, 70-84.

     Tantrum, J., Murua, A. and Stuetzle, W. (2003) Assessment and 
Pruning of Hierarchical Model Based Clustering, \emph{Proceedings of the 
ninth ACM SIGKDD international conference on Knowledge discovery and 
data mining}, Washington, D.C., 197-205.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  require(diptest)
  x <- cbind(runif(100),runif(100))
  partition <- 1+(x[,1]<0.5)
  d1 <- diptest.multi(x,partition)
  d2 <- diptest.multi(x,partition,pvalue="tantrum",M=10)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



