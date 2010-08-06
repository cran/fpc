\name{calinhara}
\alias{calinhara}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calinski-Harabasz index}
\description{
  Calinski-Harabasz index for estimating the number of clusters,
  based on an observations/variables-matrix here. A distance based
  version is available through \code{cluster.stats}.
}
\usage{
  calinhara(x,clustering,cn=max(clustering))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame.}
  \item{clustering}{vector of integers. Clustering.}
  \item{cn}{integer. Number of clusters.}
}

\value{
  Calinski-Harabasz statistic, which is 
  \code{(n-cn)*sum(diag(B))/((cn-1)*sum(diag(W)))}. B being the
  between-cluster means, 
  and W being the within-clusters covariance matrix.
}

\references{
  Calinski, R. B., and Harabasz, J. (1974) A Dendrite Method for Cluster 
Analysis, \emph{Communications in Statistics}, 3, 1-27.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{cluster.stats}}}

\examples{
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  km <- kmeans(iriss,3)
  calinhara(iriss,km$cluster)
}



\keyword{cluster}% __ONLY ONE__ keyword per line
