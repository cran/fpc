\name{cluster.magazine}
\alias{cluster.magazine}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Run many clustering methods on many numbers of clusters}
\description{
  Runs a user-specified set of clustering methods (CBI-functions, see
  \code{\link{kmeansCBI}} with several numbers of clusters on a dataset
  with unified output.
}
\usage{
cluster.magazine(data,G,diss = inherits(data, "dist"),
                             scaling=TRUE, clustermethod,
                             distmethod=rep(TRUE,length(clustermethod)),
                             ncinput=rep(TRUE,length(clustermethod)),
                             clustermethodpars,
                             trace=TRUE)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{data matrix or \code{dist}-object.}
  \item{G}{vector of integers. Numbers of clusters to consider.}
  \item{diss}{logical. If \code{TRUE}, the data matrix is assumed to be
    a distance/dissimilariy matrix, otherwise it's observations times
    variables.}
  \item{scaling}{either a logical or a numeric vector of length equal to
    the number of columns of \code{data}. If \code{FALSE}, data won't be
    scaled, otherwise \code{scaling} is passed on to \code{\link{scale}} as
    argument\code{scale}.} 
  \item{clustermethod}{vector of strings specifying names of
    CBI-functions (see \code{\link{kmeansCBI}}). These are the
    clustering methods to be applied.}
  \item{distmethod}{vector of logicals, of the same length as
    \code{clustermethod}. \code{TRUE} means that the clustering method
    operates on distances. If \code{diss=TRUE}, all entries have to be
    \code{TRUE}. Otherwise, if an entry is true, the corresponding
    method will be applied on \code{dist(data)}.}  
  \item{ncinput}{vector of logicals, of the same length as
    \code{clustermethod}. \code{TRUE} indicates that the corresponding
    clustering method requires the number of clusters as input and will
    not estimate the number of clusters itself.}
  \item{clustermethodpars}{list of the same length as
    \code{clustermethod}. Specifies parameters for all involved
    clustering methods. Its jth entry is passed to clustermethod number
    k. Can be an empty entry in case all defaults are used for a
    clustering method. The number of clusters does not need to be
    specified here.}
  \item{trace}{logical. If \code{TRUE}, some runtime information is
    printed.}
}

% \details{
% }
\value{
  List of lists comprising
  \item{output}{Two-dimensional list. The first list index i is the number
    of the clustering method (ordering as specified in
    \code{clustermethod}), the second list index j is the number of
    clusters. This stores the full output of clustermethod i run on
    number of clusters j.}
  \item{clustering}{Two-dimensional list. The first list index i is the number
    of the clustering method (ordering as specified in
    \code{clustermethod}), the second list index j is the number of
    clusters. This stores the clustering integer vector (i.e., the
    \code{partition}-component of the CBI-function, see
    \code{\link{kmeansCBI}}) of clustermethod i run on
    number of clusters j.}
  \item{noise}{Two-dimensional list. The first list index i is the number
    of the clustering method (ordering as specified in
    \code{clustermethod}), the second list index j is the number of
    clusters. List entries are single logicals. If \code{TRUE}, the
    clustering method estimated some noise, i.e., points not belonging
    to any cluster, which in the clustering vector are indicated by the
    highest number (number of clusters plus one in case that the number
    of clusters was fixed).}
  \item{othernc}{list of integer vectors of length 2. The first number is
    the number of the clustering method (the order is determined by
    argument \code{clustermethod}), the second number is the
    number of clusters for those methods that estimate the number of
    clusters themselves and estimate a number that is smaller than
    \code{min(G)} or larger than \code{max(G)}.} 
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
  \code{\link{clusterbenchstats}}, \code{\link{kmeansCBI}}
}

\examples{
  
  set.seed(20000)
  options(digits=3)
  face <- rFace(10,dMoNo=2,dNoEy=0,p=2)
  clustermethod=c("kmeansCBI","hclustCBI","hclustCBI")
# A clustering method can be used more than once, with different
# parameters
  clustermethodpars <- list()
  clustermethodpars[[2]] <- clustermethodpars[[3]] <- list()
  clustermethodpars[[2]]$method <- "complete"
  clustermethodpars[[3]]$method <- "average"
  cmf <-  cluster.magazine(face,G=2:3,clustermethod=clustermethod,
    distmethod=rep(FALSE,3),clustermethodpars=clustermethodpars)
  print(str(cmf))

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



