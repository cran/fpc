\name{findrep}
\alias{findrep}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Finding representatives for cluster border}
\description{
  Finds representative objects for the border of a cluster and the
  within-cluster variance as defined in the framework of the \code{\link{cdbw}}
  cluster validation index (and meant to be used in that context).
}
\usage{
findrep(x,xcen,clustering,cluster,r,p=ncol(x),n=nrow(x),
                    nc=sum(clustering==cluster))

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{matrix. Euclidean dataset.}
  \item{xcen}{mean vector of cluster.}
  \item{clustering}{vector of integers with length \code{=nrow(x)};
  indicating the cluster for each observation.}
  \item{cluster}{integer. Number of cluster to be treated.}
  \item{r}{integer. Number of representatives.}
  \item{p}{integer. Number of dimensions.}
  \item{n}{integer. Number of observations.}
  \item{nc}{integer. Number of observations in \code{cluster}.}
}

\value{
  List with components
  \item{repc}{vector of index of representatives (out of all
    observations).}
  \item{repx}{vector of index of representatives (out of only the
    observations in \code{cluster}).}
  \item{maxr}{number of representatives (this can be smaller than
    \code{r} if fewer pairwise different observations are in
    \code{cluster}.}
  \item{wvar}{estimated average within-cluster squared distance to mean.} 
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

\seealso{
  \code{\link{cdbw}}
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  options(digits=3)
  iriss <- as.matrix(iris[c(1:5,51:55,101:105),-5])
  irisc <- as.numeric(iris[c(1:5,51:55,101:105),5])
  findrep(iriss,colMeans(iriss),irisc,cluster=1,r=2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




