\name{stupidkcentroids}
\alias{stupidkcentroids}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stupid k-centroids random clustering}
\description{
  Picks k random centroids from given dataset and assigns every point to
  closest centroid. This is called stupid k-centroids in Hennig (2017).
}
\usage{
  stupidkcentroids(d,k)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{\code{dist}-object or dissimilarity matrix.}
  \item{k}{integer. Number of clusters.}
}

% \details{
% }
\value{
  The clustering vector (values 1 to \code{k}, length number of objects
  behind \code{d}), 
}
\references{
Hennig, C. (2017) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Proceedings of ASMDA 2017}, 501-520,
    \url{https://arxiv.org/abs/1703.09282}


}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{stupidknn}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  stupidkcentroids(dist(face),3) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



