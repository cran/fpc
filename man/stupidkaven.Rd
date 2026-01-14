\name{stupidkaven}
\alias{stupidkaven}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stupid average dissimilarity random clustering}
\description{
  Picks k random starting points from given dataset to initialise k
  clusters. Then, one by one, the point not yet assigned to any cluster
  with smallest average dissimilarity to the points of any already
  existing cluster is assigned to that
  cluster, until all points are assigned. This is a random versione of
  average linkage clustering, see
  Akhanli and Hennig (2020).
}
\usage{
  stupidkaven(d,k)
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
Akhanli, S. and Hennig, C. (2020) Calibrating and aggregating cluster
validity indexes for context-adapted comparison of clusterings.
\emph{Statistics and Computing}, 30, 1523-1544,
\url{https://link.springer.com/article/10.1007/s11222-020-09958-2}, \url{https://arxiv.org/abs/2002.01822}


}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{stupidkcentroids}}, \code{\link{stupidknn}}, \code{\link{stupidkfn}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  stupidkaven(dist(face),3) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



