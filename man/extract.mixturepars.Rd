\name{extract.mixturepars}
\alias{extract.mixturepars}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract parameters for certain components from mclust}
\description{
  Extracts parameter of certain mixture components from the output of
  \code{\link[mclust]{summary.mclustBIC}} and updates proportions so that
  they sum up to 1. 
}
\usage{
  extract.mixturepars(mclustsum,compnumbers,noise=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mclustsum}{output object of \code{\link[mclust]{summary.mclustBIC}}.}
  \item{compnumbers}{vector of integers. Numbers of mixture components.}
  \item{noise}{logical. Should be \code{TRUE} if a noise component was fitted by
    \code{\link[mclust]{mclustBIC}}.}
}

\value{
  Object as component \code{parameters} of
  \code{\link[mclust]{summary.mclustBIC}}-output, but for specified
  components only. (Orientation information from all components is kept.)
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  set.seed(98765)
  options(digits=2)
  require(mclust)
  iriss <- iris[sample(150,20),-5]
  irisBIC <- mclustBIC(iriss)
  siris <- summary(irisBIC,iriss)
  extract.mixturepars(siris,2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



