\name{cdbw}
\alias{cdbw}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{CDbw-index for cluster validation}
\description{
  CDbw-index for cluster validation, as defined in Halkidi and
  Vazirgiannis (2008), Halkidi et al. (2015).
}
\usage{
cdbw(x,clustering,r=10,s=seq(0.1,0.8,by=0.1),
                 clusterstdev=TRUE,trace=FALSE)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{something that can be coerced into a numerical
  matrix. Euclidean dataset.}
  \item{clustering}{vector of integers with length \code{=nrow(x)};
  indicating the cluster for each observation.}
  \item{r}{integer. Number of cluster border representatives.}
  \item{s}{numerical vector of shrinking factors (between 0 and 1).}
  \item{clusterstdev}{logical. If \code{TRUE}, the neighborhood radius
  for intra-cluster density is the within-cluster estimated squared
  distance from the mean of the cluster; otherwise it is the average of
  these over all clusters.}
  \item{trace}{logical. If \code{TRUE}, results are printed for the
  steps to compute the index.}
}

\value{
  List with components (see Halkidi and Vazirgiannis (2008), Halkidi et
  al. (2015) for details)
  \item{cdbw}{value of CDbw index (the higher the better).}
  \item{cohesion}{cohesion.}
  \item{compactness}{compactness.}
  \item{sep}{separation.} 
}


\references{
  Halkidi, M. and Vazirgiannis, M. (2008) A density-based cluster
  validity approach using multi-representatives. \emph{Pattern
    Recognition Letters} 29, 773-786.

  Halkidi, M., Vazirgiannis, M. and Hennig, C. (2015) Method-independent
    indices for cluster validation. In C. Hennig, M. Meila, F. Murtagh,
    R. Rocci (eds.) \emph{Handbook of Cluster Analysis}, CRC
    Press/Taylor \code{&} Francis, Boca Raton.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
  options(digits=3)
  iriss <- as.matrix(iris[c(1:5,51:55,101:105),-5])
  irisc <- as.numeric(iris[c(1:5,51:55,101:105),5])
  cdbw(iriss,irisc)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




