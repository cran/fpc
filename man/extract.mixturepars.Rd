\name{extract.mixturepars}
\alias{extract.mixturepars}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extract parameters for certain components from mclust}
\description{
  Extracts parameters of certain mixture components from the output of
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
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
  set.seed(98765)
  require(mclust)
  iriss <- iris[sample(150,20),-5]
  irisBIC <- mclustBIC(iriss,G=5,modelNames="VEV")
  siris <- summary(irisBIC,iriss)
  emp <- extract.mixturepars(siris,2)
  emp$pro
  round(emp$mean,digits=1)
  emp$variance$modelName
  round(emp$variance$scale,digits=2)
  
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



