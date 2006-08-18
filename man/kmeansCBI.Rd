\name{kmeansCBI}
\alias{kmeansCBI}
\alias{hclustCBI}
\alias{hclusttreeCBI}
\alias{disthclustCBI}
\alias{noisemclustCBI}
\alias{distnoisemclustCBI}
\alias{claraCBI}
\alias{pamkCBI}
\alias{trimkmeansCBI}
\alias{disttrimkmeansCBI}
\alias{dbscanCBI}
\alias{mahalCBI}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Interface functions for clustering methods}
\description{
  These functions provide an interface to several clustering methods
  implemented in R, for use together with the cluster stability
  assessment in \code{\link{clusterboot}} (as parameter
  \code{clustermethod}; "CBI" stands for "clusterboot interface").
  In some situations it could make sense to use them to compute a
  clustering even if you don't want to run \code{clusterboot}, because
  some of the functions contain some additional features (e.g., normal
  mixture model based clustering of dissimilarity matrices projected
  into the Euclidean space by MDS or partitioning around medoids with
  estimated number of clusters, noise/outlier identification in
  hierarchical clustering).   
}
\usage{
kmeansCBI(data,k,scaling=TRUE,runs=1,...)

hclustCBI(data,k,cut="level",method,scaling=TRUE,noisecut=0,...)

hclusttreeCBI(data,minlevel=2,method,scaling=TRUE,...)

disthclustCBI(dmatrix,k,cut="level",method,noisecut=0,...)

noisemclustCBI(data,G,emModelNames,nnk,hcmodel=NULL,Vinv=NULL)

distnoisemclustCBI(dmatrix,G,emModelNames,nnk,
                        hcmodel=NULL,Vinv=NULL,mdsmethod="classical",
                        mdsdim=4)

claraCBI(data,k,usepam=TRUE,diss=FALSE,...)

pamkCBI(data,krange=2:10,scaling=TRUE,diss=FALSE,...)

trimkmeansCBI(data,k,scaling=TRUE,trim=0.1,...)

disttrimkmeansCBI(dmatrix,k,scaling=TRUE,trim=0.1,
                         mdsmethod="classical",
                         mdsdim=4,...)

dbscanCBI(data,eps,MinPts,diss=FALSE,...)

mahalCBI(data,clustercut=0.5,...)
}
\arguments{
  \item{data}{a numeric matrix. The data
    matrix - usually a cases*variables-data matrix. \code{claraCBI},
    \code{pamkCBI} and \code{dbscanCBI} work with an
    \code{n*n}-dissimilarity matrix as well, see parameter \code{diss}.}
  \item{dmatrix}{a squared numerical dissimilarity matrix or a
    \code{dist}-object.}
  \item{k}{numeric, usually integer. In most cases, this is the number
    of clusters for methods where this is fixed. For \code{hclustCBI}
    and \code{disthclustCBI} see parameter \code{cut} below.}
  \item{scaling}{either a logical value or a numeric vector of length
    equal to the number of variables. If \code{scaling} is a numeric
    vector with length equal to the number of variables, then each
    variable is divided by the corresponding value from \code{scaling}.
    If \code{scaling} is \code{TRUE} then scaling is done by dividing
    the (centered) variables by their root-mean-square, and if
    \code{scaling} is \code{FALSE}, no scaling is done before execution.}
  \item{runs}{integer. Number of random initializations from which the
    k-means algorithm is started.}
  \item{cut}{either "level" or "number". This determines how
    \code{cutree} is used to obtain a partition from a hierarchy
    tree. \code{cut="level"} means that the tree is cut at a particular
    dissimilarity level, \code{cut="number"} means that the tree is cut
    in order to obtain a fixed number of clusters. The parameter
    \code{k} specifies the number of clusters or the dissimilarity
    level, depending on \code{cut}.}
  \item{method}{method for hierarchical clustering, see the
    documentation of \code{\link{hclust}}.}
  \item{noisecut}{numeric. All clusters of size \code{<=noisecut} in the
    \code{disthclustCBI}/\code{hclustCBI}-partition are joined and declared as
    noise/outliers.}
  \item{minlevel}{integer. \code{minlevel=1} means that all clusters in
    the tree are given out by \code{hclusttreeCBI}, including one-point
    clusters (but excluding the cluster with all
    points). \code{minlevel=2} excludes the one-point clusters.
    \code{minlevel=3} excludes the two-point cluster which has been
    merged first, and increasing the value of \code{minlevel} by 1 in
    all further steps means that the remaining earliest formed cluster
    is excluded.}
  \item{G}{vector of integers. Number of clusters or numbers of clusters
    used by
    \code{\link[mclust]{EMclust}}/\code{\link[mclust]{EMclustN}}. If
    \code{G} has more than one entry, the number of clusters is
    estimated by the BIC.}
  \item{emModelNames}{vector of string. Models for covariance matrices,
    see documentation of
    \code{\link[mclust]{EMclust}}/\code{\link[mclust]{EMclustN}}.}
  \item{nnk}{integer. Tuning constant for
    \code{\link[prabclus]{NNclean}}, which is used to estimate the
    initial noise for \code{noisemclustCBI} and
    \code{distnoisemclustCBI}. See parameter \code{k} in the
    documentation of \code{\link[prabclus]{NNclean}}. \code{nnk=0} means
    that no noise component is fitted.}
  \item{hcmodel}{string or \code{NULL}. Determines the initialization of
    the EM-algorithm for
    \code{\link[mclust]{EMclust}}/\code{\link[mclust]{EMclustN}}.
    Documented in \code{\link[mclust]{hc}}.}
  \item{Vinv}{numeric. See documentation of
    \code{\link[mclust]{EMclustN}}.}
  \item{mdsmethod}{"classical", "kruskal" or "sammon". Determines the
    multidimensional scaling method to compute Euclidean data from a
    dissimilarity matrix. See \code{\link{cmdscale}},
    \code{\link{isoMDS}} and \code{\link{sammon}}.}
  \item{mdsdim}{integer. Dimensionality of MDS solution.}
  \item{usepam}{logical. If \code{TRUE}, the function
    \code{\link[cluster]{pam}} is used for clustering, otherwise
    \code{\link[cluster]{clara}}. \code{\link{pam}} is better,
    \code{\link[cluster]{clara}} is faster.}
  \item{diss}{logical. If \code{TRUE}, \code{data} will be considered as
    a dissimilarity matrix. in \code{claraCBI}, this requires
    \code{usepam=TRUE}.}
  \item{krange}{vector of integers. Numbers of clusters to be compared.}
  \item{trim}{numeric between 0 and 1. Proportion of data points
    trimmed, i.e., assigned to noise. See
    \code{\link[trimcluster]{trimkmeans}}.}
  \item{eps}{numeric. The radius of the neighborhoods to be considered
    by \code{\link{dbscan}}.}
  \item{MinPts}{integer. How many points have to be in a neighborhood so
    that a point is considered to be a cluster seed? See documentation
    of \code{\link{dbscan}}.}
  \item{clustercut}{numeric between 0 and 1. If \code{\link{fixmahal}}
    is used for fuzzy clustering, a crisp partition is generated and
    points with cluster membership values above \code{clustercut} are
    considered as members of the corresponding cluster.}
  \item{...}{further parameters to be transferred to the original
    clustering functions (not required).}
}
\details{
  All these functions call clustering methods implemented in R to
  cluster data and to provide output in the format required by
  \code{\link{clusterboot}}. Here is a brief overview:
  \itemize{
    \item{kmeansCBI}{an interface to the function
      \code{\link{kmeansruns}} calling \code{\link{kmeans}}
      for k-means clustering. (\code{\link{kmeansruns}} allows the
      specification of several random initializations of the
      k-means algorithm.)}
    \item{hclustCBI}{an interface to the function
	\code{\link{hclust}} for agglomerative hierarchical clustering with
	noise component (see parameter \code{noisecut} above). This
	function produces a partition and assumes a cases*variables
	matrix as input.}
    \item{hclusttreeCBI}{an interface to the function
	\code{hclust} for agglomerative hierarchical clustering. This
	function gives out all clusters belonging to the hierarchy
	(upward from a certain level, see parameter \code{minlevel}
	above).}
    \item{disthclustCBI}{an interface to the function
	\code{hclust} for agglomerative hierarchical clustering with
	noise component (see parameter \code{noisecut} above). This
	function produces a partition and assumes a dissimilarity
	matrix as input.}
    \item{noisemclustCBI}{an interface to the functions
	\code{\link[mclust]{EMclust}} and
	\code{\link[mclust]{EMclustN}}, for normal mixture model based
	clustering. Warning: \code{\link[mclust]{EMclust}} and
	\code{\link[mclust]{EMclustN}} often have problems with multiple
        points. In \code{\link{clusterboot}}, it is recommended to use
	this only together with \code{multipleboot=FALSE}.}
    \item{distnoisemclustCBI}{an interface to the functions
	\code{\link[mclust]{EMclust}} and
	\code{\link[mclust]{EMclustN}}, for normal mixture model based
	clustering. This assumes a dissimilarity matrix as input and
	generates a data matrix by multidimensional scaling first.
	Warning: \code{\link[mclust]{EMclust}} and
	\code{\link[mclust]{EMclustN}} often have problems with multiple
        points. In \code{\link{clusterboot}}, it is recommended to use
	this only together with \code{multipleboot=FALSE}.}
    \item{claraCBI}{an interface to the functions
	\code{\link[cluster]{pam}} and \code{\link[cluster]{clara}}
	for partitioning around medoids.}
    \item{pamkCBI}{an interface to the function
      \code{\link{pamk}} calling \code{\link[cluster]{pam}} for
      partitioning around medoids. The number
      of clusters is estimated by the average silhouette width.}
    \item{trimkmeansCBI}{an interface to the function
	\code{\link[trimcluster]{trimkmeans}} for trimmed k-means
	clustering. This assumes a cases*variables matrix as input.}
    \item{disttrimkmeansCBI}{an interface to the function
	\code{\link[trimcluster]{trimkmeans}} for trimmed k-means
	clustering. This assumes a dissimilarity matrix as input and
	generates a data matrix by multidimensional scaling first.}
    \item{dbscanCBI}{an interface to the function
	\code{\link{dbscan}} for density based 
	clustering.}
    \item{mahalCBI}{an interface to the function
	\code{\link{fixmahal}} for fixed point
	clustering.}
  }
}

\value{
  All interface functions return a list with the following components:
  \item{result}{clustering result, usually a list with the full
	output of the clustering method (the precise format doesn't
	matter); whatever you want to use later.}
  \item{nc}{number of clusters. If some points don't belong to any
	cluster but are declared as "noise", \code{nc} includes the
	noise component, and there should be another component
	\code{nccl}, being the number of clusters not including the
	noise component.}
  \item{clusterlist}{this is a list consisting of a logical vectors
	of length of the number of data points (\code{n}) for each cluster,
	indicating whether a point is a member of this cluster
	(\code{TRUE}) or not. If a noise component is included, it
	should always be the last vector in this list.}
  \item{partition}{an integer vector of length \code{n},
	partitioning the data. If the method produces a partition, it
	should be the clustering. This component is only used for plots,
	so you could do something like \code{rep(1,n)} for
	non-partitioning methods.}
  \item{clustermethod}{a string indicating the clustering method.}      
  The output of some of the functions has further components:
  \item{nccl}{see \code{nc} above.}
  \item{nnk}{by \code{noisemclustCBI} and \code{distnoisemclustCBI},
    see above.}
  \item{initnoise}{logical vector, indicating initially estimated noise by
    \code{\link[prabclus]{NNclean}}, called by \code{noisemclustCBI}
    and \code{distnoisemclustCBI}.}
  \item{noise}{logical. \code{TRUE} if points were classified as
    noise/outliers by \code{disthclustCBI}.}
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\seealso{
  \code{\link{clusterboot}}, \code{\link{dist}},
  \code{\link{kmeans}}, \code{\link{kmeansruns}}, \code{\link{hclust}},
  \code{\link[mclust]{EMclust}}, \code{\link[mclust]{EMclustN}},
  \code{\link[cluster]{pam}},  \code{\link{pamk}},
  \code{\link[cluster]{clara}},
  \code{\link[trimcluster]{trimkmeans}}, \code{\link{dbscan}},
  \code{\link{fixmahal}}
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  dbs <- dbscanCBI(face,eps=1.5,MinPts=4)
  dhc <- disthclustCBI(dist(face),method="average",k=1.5,noisecut=2)
  table(dbs$partition,dhc$partition)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



