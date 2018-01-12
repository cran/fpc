\name{cvnn}
\alias{cvnn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cluster validation based on nearest neighbours}
\description{
  Cluster validity index based on nearest neighbours as defined in Liu
  et al. (2013) with a correction explained in Halkidi et al. (2015).
}
\usage{
cvnn(d=NULL,clusterings,k=5)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{dissimilarity matrix or \code{dist}-object.}
  \item{clusterings}{list of vectors of integers with length \code{=nrow(d)};
  indicating the cluster for each observation for several clusterings
  (list elements) to be compared.}
  \item{k}{integer. Number of nearest neighbours.}
}

\value{
  List with components (see Liu et al. (2013), Halkidi et al. (2015) for
  details)
  \item{cvnnindex}{vector of index values for the various clusterings,
    see Liu et al. (2013), the lower the better.}
  \item{sep}{vector of separation values.}
  \item{comp}{vector of compactness values.}
}

\references{
  Halkidi, M., Vazirgiannis, M. and Hennig, C. (2015) Method-independent
    indices for cluster validation. In C. Hennig, M. Meila, F. Murtagh,
    R. Rocci (eds.) \emph{Handbook of Cluster Analysis}, CRC
    Press/Taylor \code{&} Francis, Boca Raton.

  Liu, Y, Li, Z., Xiong, H., Gao, X., Wu, J. and Wu, S. (2013)
  Understanding and enhancement of internal clustering validation
  measures. \emph{IEEE Transactions on Cybernetics} 43, 982-994.
  
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  options(digits=3)
  iriss <- as.matrix(iris[c(1:10,51:55,101:105),-5])
  irisc <- as.numeric(iris[c(1:10,51:55,101:105),5])
  print(cvnn(dist(iriss),list(irisc,rep(1:4,5))))
}
\keyword{cluster}% at least one, from doc/KEYWORDS




