\name{kmeansCBI}
\alias{kmeansCBI}
\alias{hclustCBI}
\alias{hclusttreeCBI}
\alias{disthclustCBI}
\alias{disthclusttreeCBI}
\alias{noisemclustCBI}
\alias{distnoisemclustCBI}
\alias{claraCBI}
\alias{pamkCBI}
%\alias{trimkmeansCBI}
%\alias{disttrimkmeansCBI}
\alias{dbscanCBI}
\alias{mahalCBI}
\alias{mergenormCBI}
\alias{speccCBI}
\alias{tclustCBI}
\alias{pdfclustCBI}
\alias{emskewCBI}
\alias{stupidkcentroidsCBI}
\alias{stupidknnCBI}
\alias{stupidkfnCBI}
\alias{stupidkavenCBI}

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
kmeansCBI(data,krange,k,scaling=FALSE,runs=1,criterion="ch",...)

hclustCBI(data,k,cut="number",method,scaling=TRUE,noisecut=0,...)

hclusttreeCBI(data,minlevel=2,method,scaling=TRUE,...)

disthclustCBI(dmatrix,k,cut="number",method,noisecut=0,...)

% disthclusttreeCBI(dmatrix,minlevel=2,method,...)

noisemclustCBI(data,G,k,modelNames,nnk,hcmodel=NULL,Vinv=NULL,
                        summary.out=FALSE,...)

distnoisemclustCBI(dmatrix,G,k,modelNames,nnk,
                        hcmodel=NULL,Vinv=NULL,mdsmethod="classical",
                        mdsdim=4, summary.out=FALSE, points.out=FALSE,...)

claraCBI(data,k,usepam=TRUE,diss=inherits(data,"dist"),...)

pamkCBI(data,krange=2:10,k=NULL,criterion="asw", usepam=TRUE,
        scaling=FALSE,diss=inherits(data,"dist"),...)

tclustCBI(data,k,trim=0.05,...)

dbscanCBI(data,eps,MinPts,diss=inherits(data,"dist"),...)

mahalCBI(data,clustercut=0.5,...)

mergenormCBI(data, G=NULL, k=NULL, modelNames=NULL, nnk=0,
                         hcmodel = NULL,
                         Vinv = NULL, mergemethod="bhat",
                         cutoff=0.1,...)

speccCBI(data,k,...)

pdfclustCBI(data,...)

% emskewCBI(data,k,distr="mst",repeats=100,...)

stupidkcentroidsCBI(dmatrix,k,distances=TRUE)

stupidknnCBI(dmatrix,k)

stupidkfnCBI(dmatrix,k)

stupidkavenCBI(dmatrix,k)

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
    and \code{disthclustCBI} see parameter \code{cut} below. Some
    methods have a \code{k} parameter on top of a \code{G} or
    \code{krange} parameter for compatibility; \code{k} in these cases
    does not have to be specified but if it is, it is always a single
    number of clusters and overwrites \code{G} and
    \code{krange}.}
  \item{scaling}{either a logical value or a numeric vector of length
    equal to the number of variables. If \code{scaling} is a numeric
    vector with length equal to the number of variables, then each
    variable is divided by the corresponding value from \code{scaling}.
    If \code{scaling} is \code{TRUE} then scaling is done by dividing
    the (centered) variables by their root-mean-square, and if
    \code{scaling} is \code{FALSE}, no scaling is done before execution.}
  \item{runs}{integer. Number of random initializations from which the
    k-means algorithm is started.}
  \item{criterion}{\code{"ch"} or \code{"asw"}. Decides whether number
    of clusters is estimated by the Calinski-Harabasz criterion or by the
    average silhouette width.}
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
    the tree are given out by \code{hclusttreeCBI} or
    \code{disthclusttreeCBI}, including one-point
    clusters (but excluding the cluster with all
    points). \code{minlevel=2} excludes the one-point clusters.
    \code{minlevel=3} excludes the two-point cluster which has been
    merged first, and increasing the value of \code{minlevel} by 1 in
    all further steps means that the remaining earliest formed cluster
    is excluded.}
  \item{G}{vector of integers. Number of clusters or numbers of clusters
    used by
    \code{\link[mclust]{mclustBIC}}. If
    \code{G} has more than one entry, the number of clusters is
    estimated by the BIC.}
  \item{modelNames}{vector of string. Models for covariance matrices,
    see documentation of
    \code{\link[mclust]{mclustBIC}}.}
  \item{nnk}{integer. Tuning constant for
    \code{\link[prabclus]{NNclean}}, which is used to estimate the
    initial noise for \code{noisemclustCBI} and
    \code{distnoisemclustCBI}. See parameter \code{k} in the
    documentation of \code{\link[prabclus]{NNclean}}. \code{nnk=0} means
    that no noise component is fitted.}
  \item{hcmodel}{string or \code{NULL}. Determines the initialization of
    the EM-algorithm for
    \code{\link[mclust]{mclustBIC}}.
    Documented in \code{\link[mclust]{hc}}.}
  \item{Vinv}{numeric. See documentation of
    \code{\link[mclust]{mclustBIC}}.}
  \item{summary.out}{logical. If \code{TRUE}, the result of
    \code{\link[mclust]{summary.mclustBIC}} is added as component
    \code{mclustsummary} to the output of \code{noisemclustCBI} and
    \code{distnoisemclustCBI}.}
  \item{mdsmethod}{"classical", "kruskal" or "sammon". Determines the
    multidimensional scaling method to compute Euclidean data from a
    dissimilarity matrix. See \code{\link{cmdscale}},
    \code{\link[MASS]{isoMDS}} and \code{\link[MASS]{sammon}}.}
  \item{mdsdim}{integer. Dimensionality of MDS solution.}
  \item{points.out}{logical. If \code{TRUE}, the matrix of MDS points
    is added as component
    \code{points} to the output of \code{noisemclustCBI}.}
  \item{usepam}{logical. If \code{TRUE}, the function
    \code{\link[cluster]{pam}} is used for clustering, otherwise
    \code{\link[cluster]{clara}}. \code{\link[cluster]{pam}} is better,
    \code{\link[cluster]{clara}} is faster.}
  \item{diss}{logical. If \code{TRUE}, \code{data} will be considered as
    a dissimilarity matrix. In \code{claraCBI}, this requires
    \code{usepam=TRUE}.}
  \item{krange}{vector of integers. Numbers of clusters to be compared.}
  \item{trim}{numeric between 0 and 1. Proportion of data points
    trimmed, i.e., assigned to noise. See \code{tclust} in the tclust package.}
  \item{eps}{numeric. The radius of the neighborhoods to be considered
    by \code{\link{dbscan}}.}
  \item{MinPts}{integer. How many points have to be in a neighborhood so
    that a point is considered to be a cluster seed? See documentation
    of \code{\link{dbscan}}.}
  \item{clustercut}{numeric between 0 and 1. If \code{\link{fixmahal}}
    is used for fuzzy clustering, a crisp partition is generated and
    points with cluster membership values above \code{clustercut} are
    considered as members of the corresponding cluster.}
  \item{mergemethod}{method for merging Gaussians, passed on as
    \code{method} to \code{\link{mergenormals}}.}
  \item{cutoff}{numeric between 0 and 1, tuning constant for
    \code{\link{mergenormals}}.}
%  \item{distr}{one of \code{"mvn", "mvt", "msn", "mst"}. Defines the
%    family of mixtures (multivariate normal, multivariate t,
%    multivariate skew normal, or multivariate skew t). See
%    \code{\link[EMMIXskew]{EmSkew}}.}, 
%  \item{repeats}{integer. In case that \code{\link[EMMIXskew]{EmSkew}}
%    doesn't give a solution, how often should execution be repeated with
%    new random initialisations? (\code{EmSkew}'s own \code{nrandom}
%    doesn't help if an initialisation leads to a \code{NULL} output.)}
  \item{distances}{logical (only for \code{stupidkcentroidsCBI}). If
    \code{FALSE}, \code{dmatrix} is
    interpreted as cases&variables data matrix.}
  \item{...}{further parameters to be transferred to the original
    clustering functions (not required).}
}
\details{
  All these functions call clustering methods implemented in R to
  cluster data and to provide output in the format required by
  \code{\link{clusterboot}}. Here is a brief overview. For further
  details see the help pages of the involved clustering methods.
  \describe{
    \item{kmeansCBI}{an interface to the function
      \code{\link{kmeansruns}} calling \code{\link{kmeans}}
      for k-means clustering. (\code{\link{kmeansruns}} allows the
      specification of several random initializations of the
      k-means algorithm and estimation of k by the Calinski-Harabasz
      index or the average silhouette width.)}
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
%     \item{disthclusttreeCBI}{an interface to the function
% 	\code{hclust} for agglomerative hierarchical clustering. This
% 	function gives out all clusters belonging to the hierarchy
% 	(upward from a certain level, see parameter \code{minlevel}
% 	above), and assumes a dissimilarity matrix as input.}
    \item{noisemclustCBI}{an interface to the function
	\code{\link[mclust]{mclustBIC}}, for normal mixture model based
	clustering. Warning: \code{\link[mclust]{mclustBIC}} often
	has problems with multiple
        points. In \code{\link{clusterboot}}, it is recommended to use
	this together with \code{multipleboot=FALSE}.}	
    \item{distnoisemclustCBI}{an interface to the function
	\code{\link[mclust]{mclustBIC}} for normal mixture model based
	clustering. This assumes a dissimilarity matrix as input and
	generates a data matrix by multidimensional scaling first.
	Warning: \code{\link[mclust]{mclustBIC}} often has
	problems with multiple
        points. In \code{\link{clusterboot}}, it is recommended to use
	this together with \code{multipleboot=FALSE}.}
    \item{claraCBI}{an interface to the functions
	\code{\link[cluster]{pam}} and \code{\link[cluster]{clara}}
	for partitioning around medoids.}
    \item{pamkCBI}{an interface to the function
      \code{\link{pamk}} calling \code{\link[cluster]{pam}} for
      partitioning around medoids. The number
      of clusters is estimated by the Calinski-Harabasz index or by the
      average silhouette width.}
    \item{tclustCBI}{an interface to the function
	\code{tclust} in the tclust package for trimmed Gaussian 
	clustering. This assumes a cases*variables matrix as input.}
%	
%        NOTE: This package is currently only available in CRAN as
%	archived version. Therefore I cannot currently offer the
%	\code{tclustCBI}-function in \code{fpc}. The code for the
%	function is below in the Examples-Section, so if you need it,
%	run that code first.}
%    \item{disttrimkmeansCBI}{an interface to the function
%	\code{\link[trimcluster]{trimkmeans}} for trimmed k-means
%	clustering. This assumes a dissimilarity matrix as input and
%	generates a data matrix by multidimensional scaling first.}
    \item{dbscanCBI}{an interface to the function
	\code{\link{dbscan}} for density based 
	clustering.}
    \item{mahalCBI}{an interface to the function
	\code{\link{fixmahal}} for fixed point
	clustering. This assumes a cases*variables matrix as input.}
    \item{mergenormCBI}{an interface to the function
      \code{\link{mergenormals}} for clustering by merging Gaussian
      mixture components. Unlike \code{\link{mergenormals}}, \code{mergenormCBI}
      includes the computation of the initial Gaussian mixture.
      This assumes a cases*variables matrix as input.
    }
    \item{speccCBI}{an interface to the function
      \code{\link[kernlab]{specc}} for spectral clustering. See
      the \code{\link[kernlab]{specc}} help page for additional tuning
      parameters. This assumes a cases*variables matrix as input.}
    \item{pdfclustCBI}{an interface to the function
      \code{\link[pdfCluster]{pdfCluster}} for density-based clustering. See
      the \code{\link[pdfCluster]{pdfCluster}} help page for additional tuning
      parameters. This assumes a cases*variables matrix as input.}
%    \item{emskewCBI}{an interface to the function
%      \code{\link[EMMIXskew]{EmSkew}} for clustering with the
%      EM-algorithm based on Gaussian, skew Gaussian, t or skew-t
%      mixtures. See
%      help page of \code{\link[EMMIXskew]{EmSkew}}. This assumes a
%      cases*variables matrix as input. Note that by September 2020,
%      package \code{EMMIXskew} is not available on CRAN but only
%      in the CRAN archives; CRAN states that it needs an update.} 
    \item{stupidkcentroidsCBI}{an interface to the function
      \code{stupidkcentroids} for random centroid-based clustering. See
      the \code{\link{stupidkcentroids}} help page. This can have a
      distance matrix as well as a cases*variables matrix as input, see
      parameter \code{distances}.}
    \item{stupidknnCBI}{an interface to the function
      \code{stupidknn} for random nearest neighbour clustering. See
      the \code{\link{stupidknn}} help page. This assumes a
      distance matrix as input.}
    \item{stupidkfnCBI}{an interface to the function
      \code{stupidkfn} for random farthest neighbour clustering. See
      the \code{\link{stupidkfn}} help page. This assumes a
      distance matrix as input.}
    \item{stupidkavenCBI}{an interface to the function
      \code{stupidkaven} for random average dissimilarity clustering. See
      the \code{\link{stupidkaven}} help page. This assumes a
      distance matrix as input.}
  }
}

\value{
  All interface functions return a list with the following components
  (there may be some more, see \code{summary.out} and \code{points.out}
  above):
  \item{result}{clustering result, usually a list with the full
	output of the clustering method (the precise format doesn't
	matter); whatever you want to use later.}
  \item{nc}{number of clusters. If some points don't belong to any
	cluster, these are declared "noise". \code{nc} includes the
	"noise cluster", and there should be another component
	\code{nccl}, being the number of clusters not including the
	noise cluster.}
  \item{clusterlist}{this is a list consisting of a logical vectors
	of length of the number of data points (\code{n}) for each cluster,
	indicating whether a point is a member of this cluster
	(\code{TRUE}) or not. If a noise cluster is included, it
	should always be the last vector in this list.}
  \item{partition}{an integer vector of length \code{n},
	partitioning the data. If the method produces a partition, it
	should be the clustering. This component is only used for plots,
	so you could do something like \code{rep(1,n)} for
	non-partitioning methods. If a noise cluster is included,
        \code{nc=nccl+1} and the noise cluster is cluster no. \code{nc}.}
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
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{clusterboot}}, \code{\link{dist}},
  \code{\link{kmeans}}, \code{\link{kmeansruns}}, \code{\link{hclust}},
  \code{\link[mclust]{mclustBIC}}, 
  \code{\link[cluster]{pam}},  \code{\link{pamk}},
  \code{\link[cluster]{clara}},
  \code{\link{dbscan}},
  \code{\link{fixmahal}},
  \code{\link[tclust]{tclust}}, \code{\link[pdfCluster]{pdfCluster}}
%  \code{\link[EMMIXskew]{EmSkew}}
}
\examples{
  options(digits=3)
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  dbs <- dbscanCBI(face,eps=1.5,MinPts=4)
  dhc <- disthclustCBI(dist(face),method="average",k=1.5,noisecut=2)
  table(dbs$partition,dhc$partition)
  dm <- mergenormCBI(face,G=10,modelNames="EEE",nnk=2)
  dtc <- tclustCBI(face,6,trim=0.1,restr.fact=500)
  table(dm$partition,dtc$partition)

}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}