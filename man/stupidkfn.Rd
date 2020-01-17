\name{stupidkfn}
\alias{stupidkfn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stupid farthest neighbour random clustering}
\description{
  Picks k random starting points from given dataset to initialise k
  clusters. Then, one by one, a point not yet assigned to any cluster
  is assigned to that
  cluster, until all points are assigned. The point/cluster pair to be
  used is picked according to the smallest distance of a point to the
  farthest point to it in any of the already existing clusters as in
  complete linkage clustering, see
  Akhanli and Hennig (2020).
}
\usage{
  stupidkfn(d,k)
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
On arxiv from February 2020. 


}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}

\seealso{
  \code{\link{stupidkcentroids}}, \code{\link{stupidknn}}, \code{\link{stupidkaven}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  stupidkfn(dist(face),3) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



