\name{dudahart2}
\alias{dudahart2}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Duda-Hart test for splitting}
\description{
  Duda-Hart test for whether a data set should be split into two
  clusters. 
}
\usage{
  dudahart2(x,clustering,alpha=0.001)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame.}
  \item{clustering}{vector of integers. Clustering into two clusters.}
  \item{alpha}{numeric between 0 and 1. Significance level (recommended
    to be small if this is used for estimating the number of clusters).}
}

\value{
  A list with components
  \item{p.value}{p-value against null hypothesis of homogemeity.}
  \item{dh}{ratio of within-cluster sum of squares for two clusters and
    overall sum of squares.}
  \item{compare}{critical value for \code{dh} at level \code{alpha}.}
  \item{cluster1}{\code{FALSE} if the null hypothesis of homogemeity is
    rejected.}
  \item{alpha}{see above.}
  \item{z}{\code{1-alpha}-quantile of a standard Gaussian.}
}

\references{
  Duda, R. O. and Hart, P. E. (1973) \emph{Pattern Classification and
  Scene Analysis}. Wiley, New York.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{cluster.stats}}}

\examples{
  options(digits=2)
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  km <- kmeans(iriss,2)
  dudahart2(iriss,km$cluster)
}

\keyword{cluster}% __ONLY ONE__ keyword per line
