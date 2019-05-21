\name{cgrestandard}
\alias{cgrestandard}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Standardise cluster validation statistics by random clustering results}
\description{
  Standardises cluster validity statistics as produced by
  \code{\link{clustatsum}} relative to results that were achieved by
  random clusterings on the same data by
  \code{\link{randomclustersim}}. The aim is to make differences between
  values comparable between indexes, see Hennig (2017).

  This is mainly for use within \code{\link{clusterbenchstats}}.
}
\usage{
cgrestandard(clusum,clusim,G,percentage=FALSE,
                               useallmethods=FALSE,
                             useallg=FALSE, othernc=list())
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{clusum}{object of class "valstat", see \code{\link{clusterbenchstats}}.}
  \item{clusim}{list; output object of \code{\link{randomclustersim}},
    see there.} 
  \item{G}{vector of integers. Numbers of clusters to consider.}
  \item{percentage}{logical. If \code{FALSE}, standardisation is done to
    mean zero and standard deviation 1 using the random clusterings. If
    \code{TRUE}, the output is the percentage of simulated values below
    the result (more precisely, this number plus one divided by the
    total plus one).}
  \item{useallmethods}{logical. If \code{FALSE}, only random clustering
    results from \code{clusim} are used for standardisation. If
    \code{TRUE}, also clustering results from other methods as given in
    \code{clusum} are used.}
  \item{useallg}{logical. If \code{TRUE}, standardisation uses results
    from all numbers of clusters in \code{G}. If \code{FALSE},
    standardisation of results for a specific number of cluster only
    uses results from that number of clusters.} 
  \item{othernc}{list of integer vectors of length 2. This allows the
    incorporation of methods that bring forth other numbers of clusters
    than those in \code{G}, for example because a method may have
    automatically estimated a number of clusters. The first number is
    the number of the clustering method (the order is determined by
    argument \code{clustermethod} in
    \code{\link{clusterbenchstats}}), the second number is the
    number of clusters. Results specified here are only standardised in
    \code{useallg=TRUE}.}
}

\details{
  \code{cgrestandard} will add a statistic named \code{dmode} to the
  input set of validation statistics, which is defined as
  \code{0.75*dindex+0.25*highdgap}, aggregating these two closely
  related statistics, see \code{\link{clustatsum}}.
}

\value{
  List of class \code{"valstat"}, see
  \code{\link{valstat.object}}, with standardised results as
  explained above.
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
  \code{\link{valstat.object}}, \code{\link{clusterbenchstats}}, \code{\link{stupidkcentroids}}, \code{\link{stupidknn}}, \code{\link{clustatsum}}
}

\examples{
  
  set.seed(20000)
  options(digits=3)
  face <- rFace(10,dMoNo=2,dNoEy=0,p=2)
  dif <- dist(face)
  clusum <- list()
  clusum[[2]] <- list()
  cl12 <- kmeansCBI(face,2)
  cl13 <- kmeansCBI(face,3)
  cl22 <- claraCBI(face,2)
  cl23 <- claraCBI(face,2)
  ccl12 <- clustatsum(dif,cl12$partition)
  ccl13 <- clustatsum(dif,cl13$partition)
  ccl22 <- clustatsum(dif,cl22$partition)
  ccl23 <- clustatsum(dif,cl23$partition)
  clusum[[1]] <- list()
  clusum[[1]][[2]] <- ccl12
  clusum[[1]][[3]] <- ccl13
  clusum[[2]][[2]] <- ccl22
  clusum[[2]][[3]] <- ccl23
  clusum$maxG <- 3
  clusum$minG <- 2
  clusum$method <- c("kmeansCBI","claraCBI")
  clusum$name <- c("kmeansCBI","claraCBI")
  clusim <- randomclustersim(dist(face),G=2:3,nnruns=2,kmruns=2,monitor=FALSE)
  cgr <- cgrestandard(clusum,clusim,2:3)
  cgr2 <- cgrestandard(clusum,clusim,2:3,useallg=TRUE)
  cgr3 <- cgrestandard(clusum,clusim,2:3,percentage=TRUE)
  print(str(cgr))
  print(str(cgr2))
  print(str(cgr3))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



