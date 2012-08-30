\name{dipp.tantrum}
\alias{dipp.tantrum}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Simulates p-value for dip test}
\description{
  Simulates p-value for dip test (see \code{\link[diptest]{dip}})
  in the way suggested by Tantrum, Murua and Stuetzle (2003) from the
  clostest unimodal distribution determined by kernel density estimation
  with bandwith chosen so that the density just becomes unimodal. This is
  less conservative (and in fact sometimes anti-conservative) than the
  values from \code{\link[diptest]{dip.test}}.
}
\usage{
  dipp.tantrum(xdata,d,M=100)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{numeric vector. One-dimensional dataset.}
  \item{d}{numeric. Value of dip statistic.}
  \item{M}{integer. Number of artificial datasets generated in order to
    estimate the p-value.}
}

\value{
  List with components
  \item{p.value}{approximated p-value.}
  \item{bw}{borderline unimodality bandwith in \code{\link{density}}
    with default settings.}
  \item{dv}{vector of dip statistic values from simulated artificial data.}
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
# not run, requires package diptest
#  x <- runif(100)
#  d <- dip(x)
#  dt <- dipp.tantrum(x,d,M=10)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
% \keyword{multivariate}



