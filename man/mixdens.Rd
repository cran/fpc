\name{mixdens}
\alias{mixdens}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density of multivariate Gaussian mixture, mclust parameterisation}
\description{
  Computes density values for data from a mixture of multivariate Gaussian
  distributions with parameters based on the way models are specified
  and parameters are stored in package mclust.
}
\usage{
  mixdens(modelName,data,parameters)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{modelName}{an mclust model name.
    See \code{\link[mclust]{mclustModelNames}}.}
  \item{data}{data matrix; density values are computed for every
    observation (row).}
  \item{parameters}{parameters of Gaussian mixture in the format used in
    the output of \code{\link[mclust]{summary.mclustBIC}}.}
}

\value{
  Vector of density values for the observations.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  set.seed(98765)
  require(mclust)
  iriss <- iris[sample(150,20),-5]
  irisBIC <- mclustBIC(iriss)
  siris <- summary(irisBIC,iriss)
  round(mixdens(siris$modelName,iriss,siris$parameters),digits=2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



