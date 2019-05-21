\name{stupidknn}
\alias{stupidknn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stupid nearest neighbour random clustering}
\description{
  Picks k random starting points from given dataset to initialise k
  clusters. Then, one by one, the point not yet assigned to any cluster
  that is closest to an already assigned point is assigned to that
  cluster, until all points are assigned. This is called stupid nearest
  neighbour clustering in Hennig (2017).
}
\usage{
  stupidknn(d,k)
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
  \code{\link{stupidkcentroids}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  stupidknn(dist(face),3) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



