\name{stupidkcentroids}
\alias{stupidkcentroids}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stupid k-centroids random clustering}
\description{
  Picks k random centroids from given dataset and assigns every point to
  closest centroid. This is called stupid k-centroids in Hennig (2019).
}
\usage{
  stupidkcentroids(xdata, k, distances = inherits(xdata, "dist"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{cases*variables data, \code{dist}-object or dissimilarity
    matrix, see \code{distances}.}
  \item{k}{integer. Number of clusters.}
  \item{distances}{logical. If \code{TRUE}, \code{xdata} is interpreted
    as distances.}
}

% \details{
% }
\value{
  The clustering vector (values 1 to \code{k}, length number of objects
  behind \code{xdata}), 
}
\references{

Hennig, C. (2019) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Data Analysis and Applications 1: Clustering and Regression,
  Modeling-estimating, Forecasting and Data Mining, Volume 2}, Wiley,
  New York 1-24,
    \url{https://arxiv.org/abs/1703.09282}

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
  \code{\link{stupidknn}}, \code{\link{stupidkfn}}, \code{\link{stupidkaven}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(200,dMoNo=2,dNoEy=0,p=2)
  stupidkcentroids(dist(face),3) 
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



