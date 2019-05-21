\name{clusterbenchstats}
\alias{clusterbenchstats}
\alias{print.clusterbenchstats}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Run and validate many clusterings}
\description{
  This runs the methodology explained in Hennig (2017). It runs a
  user-specified set of clustering methods (CBI-functions, see
  \code{\link{kmeansCBI}}) with several numbers of clusters on a dataset,
  and computes many cluster validation indexes. In order to explore the
  variation of these indexes, random clusterings on the data are
  generated, and validation indexes are standardised by use of the
  random clusterings in order to make them comparable and differences
  between values interpretable.

  The function \code{\link{print.valstat}} can be used to provide
  weights for the cluster
  validation statistics, and will then compute a weighted validation index
  that can be used to compare all clusterings. 
}
\usage{
clusterbenchstats(data,G,diss = inherits(data, "dist"),
                                  scaling=TRUE, clustermethod,
                                  methodnames=clustermethod,
                              distmethod=rep(TRUE,length(clustermethod)),
                              ncinput=rep(TRUE,length(clustermethod)),
                              clustermethodpars,
                              npstats=FALSE,
                              trace=TRUE,
                              pamcrit=TRUE,snnk=2,
                              dnnk=2,
                              nnruns=100,kmruns=100,
                              multicore=FALSE,cores=detectCores()-1,
                              useallmethods=TRUE,
                              useallg=FALSE,...)

\method{print}{clusterbenchstats}(x,...)

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
  \item{methodnames}{vector of strings with user-chosen names for
    clustering methods, one for every method in
    \code{clustermethod}. These can be used to distinguish different methods
    run by the same CBI-function but with
    different parameter values such as complete and average linkage for
    \code{\link{hclustCBI}}.}
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
  \item{npstats}{logical. If \code{TRUE}, \code{\link{distrsimilarity}}
    is called and the two validity statistics computed there are
    added. These require \code{diss=FALSE}.}
  \item{trace}{logical. If \code{TRUE}, some runtime information is
    printed.}
  \item{pamcrit}{logical. If \code{TRUE}, the average distance of points
    to their respective cluster centroids is computed (criterion of the
    PAM clustering method, validation criterion \code{pamc}); centroids
    are chosen so that they minimise
    this criterion for the given clustering. Passed on to
    \code{\link{cqcluster.stats}}.}
  \item{snnk}{integer. Number of neighbours used in coefficient of
    variation of distance to nearest within cluster neighbour, the
    \code{cvnnd}-statistic  (clusters
    with \code{nnk} or fewer points are ignored for this). Passed on to
    \code{\link{cqcluster.stats}}.}
    \item{dnnk}{integer. Number of nearest neighbors to use for
    dissimilarity to the uniform in case that \code{npstats=TRUE};
    \code{nnk}-argument  to be passed on to \code{\link{distrsimilarity}}.}
  \item{nnruns}{integer. Number of runs of \code{\link{stupidknn}}
    (random clusterings).}
  \item{kmruns}{integer. Number of runs of
    \code{\link{stupidkcentroids}} (random clusterings).}
  \item{multicore}{logical. If \code{TRUE}, parallel computing is used
    through the function \code{\link{mclapply}} from package
    \code{parallel}; read warnings there if you intend to use this; it
    won't work on Windows.}
  \item{cores}{integer. Number of cores for parallelisation.}
  \item{useallmethods}{logical, to be passed on to
    \code{\link{cgrestandard}}. If \code{FALSE}, only random clustering
    results are used for standardisation. If
    \code{TRUE}, clustering results from all methods are used.}
  \item{useallg}{logical to be passed on to
    \code{\link{cgrestandard}}. If \code{TRUE}, standardisation uses results
    from all numbers of clusters in \code{G}. If \code{FALSE},
    standardisation of results for a specific number of cluster only
    uses results from that number of clusters.}
  \item{...}{further arguments to be passed on to
    \code{\link{cqcluster.stats}} through \code{\link{clustatsum}} (no
    effect in \code{print.clusterbenchstats}).}
  \item{x}{object of class \code{"clusterbenchstats"}.}


}

\note{
  This may require a lot of computing time and also memory for datasets
  that are not small, as most indexes require computation and storage of
  distances. 
}

\value{
  The output of \code{clusterbenchstats} is a
  big list of lists comprising lists \code{cm, stat, sim, qstat,
    sstat, statistics}
  
  \item{cm}{output object of \code{\link{cluster.magazine}}, see there
    for details. Clustering of all methods and numbers of clusters on
    the dataset \code{data}.}. 
  \item{stat}{object of class \code{"valstat"}, see
    \code{\link{valstat.object}} for details. Unstandardised cluster
    validation statistics.}
  \item{sim}{output object of \code{\link{randomclustersim}}, see there.
    validity indexes from random clusterings used for standardisation of
    validation statistics on \code{data}.}
  \item{qstat}{object of class \code{"valstat"}, see
    \code{\link{valstat.object}} for details. Cluster validation
    statistics standardised by random clusterings, output of
    \code{\link{cgrestandard}} based on percentages, i.e., with
    \code{percentage=TRUE}.}
  \item{sstat}{object of class \code{"valstat"}, see
    \code{\link{valstat.object}} for details. Cluster validation
    statistics standardised by random clusterings, output of
    \code{\link{cgrestandard}} based on mean and standard deviation,
    i.e., with \code{percentage=FALSE}.}

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
  \code{\link{valstat.object}},
  \code{\link{cluster.magazine}}, \code{\link{kmeansCBI}},
  \code{\link{cqcluster.stats}}, \code{\link{clustatsum}},
  \code{\link{cgrestandard}}  
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
  methodname <- c("kmeans","complete","average")
  cbs <-  clusterbenchstats(face,G=2:3,clustermethod=clustermethod,
    methodname=methodname,distmethod=rep(FALSE,3),
    clustermethodpars=clustermethodpars,nnruns=2,kmruns=2)
  print(cbs)
  print(cbs$qstat,aggregate=TRUE,weights=c(1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1))
# The weights are weights for the validation statistics ordered as in
# cbs$qstat$statistics for computation of an aggregated index, see
# ?print.valstat.

}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line


