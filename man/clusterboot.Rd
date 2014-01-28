\name{clusterboot}
\alias{clusterboot}
\alias{print.clboot}
\alias{plot.clboot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Clusterwise cluster stability assessment by resampling}
\description{
  Assessment of the clusterwise stability of a clustering of data, which
  can be cases*variables or dissimilarity data. The data is resampled
  using several schemes (bootstrap, subsetting, jittering, replacement
  of points by noise) and the Jaccard similarities of the original
  clusters to the most similar clusters in the resampled data are
  computed. The mean over these similarities is used as an index of the
  stability of a cluster (other statistics can be computed as well). The
  methods are described in Hennig (2007).

  \code{clusterboot} is an integrated function that computes the
  clustering as well, using interface functions for various
  clustering methods implemented in R (several interface functions are
  provided, but you can
  implement further ones for your favourite clustering method). See the
  documentation of the input parameter \code{clustermethod} below.

  Quite general clustering methods are possible, i.e. methods estimating
  or fixing the number of clusters, methods producing overlapping
  clusters or not assigning all cases to clusters (but declaring them as
  "noise"). Fuzzy clusterings cannot be processed and have to be
  transformed to crisp clusterings by the interface function.  
}
\usage{
clusterboot(data,B=100, distances=(class(data)=="dist"),
                        bootmethod="boot",
                        bscompare=TRUE, 
                        multipleboot=FALSE,
                        jittertuning=0.05, noisetuning=c(0.05,4),
                        subtuning=floor(nrow(data)/2),
                        clustermethod,noisemethod=FALSE,count=TRUE,
                        showplots=FALSE,dissolution=0.5,
                        recover=0.75,seed=NULL,...)

\method{print}{clboot}(x,statistics=c("mean","dissolution","recovery"),...)

\method{plot}{clboot}(x,xlim=c(0,1),breaks=seq(0,1,by=0.05),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{something that can be coerced into a matrix. The data
    matrix - either an \code{n*p}-data matrix (or data frame) or an
    \code{n*n}-dissimilarity matrix (or \code{dist}-object).}
  \item{B}{integer. Number of resampling runs for each scheme, see
    \code{bootmethod}.}
  \item{distances}{logical. If \code{TRUE}, the data is interpreted as
    dissimilarity matrix. If \code{data} is a \code{dist}-object,
    \code{distances=TRUE} automatically, otherwise
    \code{distances=FALSE} by default. This means that you have to set
    it to \code{TRUE} manually if \code{data} is a dissimilarity matrix.}
  \item{bootmethod}{vector of strings, defining the methods used for
    resampling. Possible methods:
   
      \code{"boot"}: nonparametric bootstrap (precise behaviour is
      controlled by parameters \code{bscompare} and
      \code{multipleboot}).
      
      \code{"subset"}: selecting random subsets from the dataset. Size
      determined by \code{subtuning}.
      
      \code{"noise"}: replacing a certain percentage of the points by
      random noise, see \code{noisetuning} (note that this will not
      work if within-cluster .
      
      \code{"jitter"} add random noise to all points, see
      \code{jittertuning}. (This didn't perform well in Hennig (2007),
      but you may want to get your own experience.)
      
      \code{"bojit"} nonparametric bootstrap first, and then adding
      noise to the points, see \code{jittertuning}.
      
    \strong{Important:} only the methods \code{"boot"} and
    \code{"subset"} work with dissimilarity data!

    The results in Hennig (2007) indicate that \code{"boot"} is
    generally informative and often quite similar to \code{"subset"} and
    \code{"bojit"}, while \code{"noise"} sometimes provides different
    information. Therefore the default (for \code{distances=FALSE}) is
    to use \code{"boot"} and \code{"noise"}. However, some clustering
    methods may have problems with multiple points, which can be solved
    by using \code{"bojit"} or \code{"subset"} instead of \code{"boot"} or by
    \code{multipleboot=FALSE} below.}
  \item{bscompare}{logical. If \code{TRUE}, multiple points in the
    bootstrap sample are taken into account to compute the Jaccard
    similarity to the original clusters (which are represented by their
    "bootstrap versions", i.e., the
    points of the original cluster which also occur in the bootstrap
    sample). If a point was drawn more than once, it is in the "bootstrap
    version" of the original cluster more than once, too, if
    \code{bscompare=TRUE}. Otherwise multiple points are
    ignored for the computation of the Jaccard similarities. If
    \code{multipleboot=FALSE}, it doesn't make a difference.} 
  \item{multipleboot}{logical. If \code{FALSE}, all points drawn more
    than once in the bootstrap draw are only used once in the bootstrap
    samples.}
  \item{jittertuning}{positive numeric. Tuning for the
    \code{"jitter"}-method. The noise distribution for
    jittering is a normal distribution with zero mean. The covariance
    matrix has the same Eigenvectors as that of the original
    data set, but the standard deviation along the principal directions is
    determined by the \code{jittertuning}-quantile of the distances
    between neighboring points projected along these directions.}
  \item{noisetuning}{A vector of two positive numerics. Tuning for the
    \code{"noise"}-method. The first component determines the
    probability that a point is replaced by noise. Noise is generated by
    a uniform distribution on a hyperrectangle along the principal
    directions of the original data set, ranging from
    \code{-noisetuning[2]} to \code{noisetuning[2]} times the standard
    deviation of the data set along the respective direction. Note that
    only points not replaced by noise are considered for the computation
    of Jaccard similarities.}
  \item{subtuning}{integer. Size of subsets for \code{"subset"}.}
  \item{clustermethod}{an interface function (the function name, not a
    string containing the name, has to be provided!). This defines the
    clustering method. See the "Details"-section for a list of available
    interface functions and guidelines how to write your own ones. 
  }
  \item{noisemethod}{logical. If \code{TRUE}, the last cluster is
    regarded as "noise component", which means that for computing the Jaccard
    similarity, it is not treated as a cluster. The noise component of
    the original clustering is only compared with the noise component of
    the clustering of the resampled data. This means that in the
    \code{clusterboot}-output (and plot), if points were assigned to the
    noise component, the last cluster number refers to it, and its
    Jaccard similarity values refer to comparisons with estimated noise
    components in resampled datasets only.
    (Some cluster methods such as
    trimmed k-means and \code{\link[mclust]{mclustBIC}} produce such noise
    components.)}
  \item{count}{logical. If \code{TRUE}, the resampling runs are counted
    on the screen.}    
  \item{showplots}{logical. If \code{TRUE}, a plot of the first two
    dimensions of the resampled data set (or the classical MDS solution
    for dissimilarity data) is shown for every resampling run. The last
    plot shows the original data set.}
  \item{dissolution}{numeric between 0 and 1. If the Jaccard similarity
    between the resampling version of the original cluster and the most
    similar cluster on the resampled data is smaller or equal to this
    value, the cluster is considered as "dissolved". Numbers of
    dissolved clusters are recorded.}
  \item{recover}{numeric between 0 and 1. If the Jaccard similarity
    between the resampling version of the original cluster and the most
    similar cluster on the resampled data is larger than  this
    value, the cluster is considered as "successfully recovered". Numbers of
    recovered clusters are recorded.}
  \item{seed}{integer. Seed for random generator (fed into
    \code{set.seed}) to make results reproducible. If \code{NULL},
    results depend on chance.}
  \item{...}{additional parameters for the clustermethods called by
    \code{clusterboot}. No effect in \code{print.clboot} and
    \code{plot.clboot}.}
  \item{x}{object of class \code{clboot}.}
  \item{statistics}{specifies in \code{print.clboot},
    which of the three clusterwise Jaccard
    similarity statistics \code{"mean"}, \code{"dissolution"} (number of
    times the cluster has been dissolved) and \code{"recovery"} (number
    of times a cluster has been successfully recovered) is printed.}
  \item{xlim}{transferred to \code{hist}.}
  \item{breaks}{transferred to \code{hist}.}
}

\details{
  Here are some guidelines for interpretation.
  There is some theoretical justification to consider a Jaccard
  similarity value smaller or equal to 0.5 as an indication of a
  "dissolved cluster", see Hennig (2008). Generally, a valid, stable
  cluster should yield a mean Jaccard similarity value of 0.75 or more.
  Between 0.6 and 0.75, clusters may be considered as indicating
  patterns in the data, but which points exactly should belong to these
  clusters is highly doubtful. Below average Jaccard values of 0.6, clusters
  should not be trusted. "Highly stable" clusters should yield average
  Jaccard similarities of 0.85 and above. All of this refers to
  bootstrap; for the other resampling schemes it depends on the tuning
  constants, though their default values should grant similar
  interpretations in most cases.

  While \code{B=100} is recommended, smaller run numbers could give
  quite informative results as well, if computation times become too high.

  Note that the stability of a cluster is assessed, but
  stability is not the only important validity criterion - clusters
  obtained by very inflexible clustering methods may be stable but not
  valid, as discussed in Hennig (2007).
  See \code{\link{plotcluster}} for graphical cluster validation.

  Information about interface functions for clustering methods:
  
    The following interface functions are currently
    implemented (in the present package; note that almost all of these
    functions require the specification of some control parameters, so
    if you use one of them, look up their common help page
    \code{\link{kmeansCBI}}) first:
    \describe{
      \item{kmeansCBI}{an interface to the function
	\code{\link{kmeans}} for k-means clustering. This assumes a
	cases*variables matrix as input.}
      \item{hclustCBI}{an interface to the function
	\code{hclust} for agglomerative hierarchical clustering with
	optional noise component. This
	function produces a partition and assumes a cases*variables
	matrix as input.}
      \item{hclusttreeCBI}{an interface to the function
	\code{hclust} for agglomerative hierarchical clustering. This
	function produces a tree (not only a partition; therefore the
	number of clusters can be huge!) and assumes a cases*variables
	matrix as input.}
      \item{disthclustCBI}{an interface to the function
	\code{hclust} for agglomerative hierarchical clustering with
	optional noise component. This
	function produces a partition and assumes a dissimilarity
	matrix as input.}
      \item{noisemclustCBI}{an interface to the function
	\code{\link[mclust]{mclustBIC}} for normal mixture model based
	clustering. This assumes a cases*variables matrix as
	input. Warning: \code{\link[mclust]{mclustBIC}} sometimes has
	problems with multiple
        points. It is recommended to use this only together with
	\code{multipleboot=FALSE}.}
      \item{distnoisemclustCBI}{an interface to the function
	\code{\link[mclust]{mclustBIC}} for normal mixture model based
	clustering. This assumes a dissimilarity matrix as input and
	generates a data matrix by multidimensional scaling first.
	Warning: \code{\link[mclust]{mclustBIC}} sometimes has
	problems with multiple
        points. It is recommended to use this only together with
	\code{multipleboot=FALSE}.}
      \item{claraCBI}{an interface to the functions
	\code{\link[cluster]{pam}} and \code{\link[cluster]{clara}}
	for partitioning around medoids. This can be used with
	cases*variables as well as dissimilarity matrices as input.}
      \item{pamkCBI}{an interface to the function
	\code{\link{pamk}} for partitioning around medoids. The number
	of cluster is estimated by the average silhouette width.
	This can be used with
	cases*variables as well as dissimilarity matrices as input.}
      \item{trimkmeansCBI}{an interface to the function
	\code{\link[trimcluster]{trimkmeans}} for trimmed k-means
	clustering. This assumes a cases*variables matrix as input.}
      \item{tclustCBI}{an interface to the function
	\code{tclust} in the tclust library for trimmed Gaussian 
	clustering. This assumes a cases*variables matrix as input. Note
        that this function is not currently provided because the tclust
        package is only available in the CRAN archives, but the code is
        in the Examples-section of the \code{\link{kmeansCBI}}-help page.}
      \item{disttrimkmeansCBI}{an interface to the function
	\code{\link[trimcluster]{trimkmeans}} for trimmed k-means
	clustering. This assumes a dissimilarity matrix as input and
	generates a data matrix by multidimensional scaling first.}
      \item{dbscanCBI}{an interface to the function
	\code{\link{dbscan}} for density based 
	clustering. This can be used with
	cases*variables as well as dissimilarity matrices as input..}
      \item{mahalCBI}{an interface to the function
	\code{\link{fixmahal}} for fixed point
	clustering. This assumes a cases*variables matrix as input.}
      \item{mergenormCBI}{an interface to the function
      \code{\link{mergenormals}} for clustering by merging Gaussian
      mixture components.}
      \item{speccCBI}{an interface to the function
      \code{\link[kernlab]{specc}} for spectral clustering.}
    }
    You can write your own interface function. The first argument of an
    interface function should always be a data matrix (of class
    "matrix", but it may be a symmetrical dissimilarity matrix). Further
    arguments can be tuning constants for the clustering method. The
    output of an interface function should be a list containing (at
    least) the following components:
    \describe{
      \item{result}{clustering result, usually a list with the full
	output of the clustering method (the precise format doesn't
	matter); whatever you want to use later.}
      \item{nc}{number of clusters. If some points don't belong to any
	cluster but are declared as "noise", \code{nc} includes the
	noise component, and there should be another component
	\code{nccl}, being the number of clusters not including the
	noise component (note that it is not mandatory to define a noise
	component if not all points are assigned to clusters, but if you
	do it, the stability of the noise component is assessed as
	well.)}
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
    }      
  
}

\value{
  \code{clusterboot} returns an object of class \code{"clboot"}, which
  is a list with components
  \code{result, partition, nc, clustermethod, B, noisemethod, bootmethod,
    multipleboot, dissolution, recover, bootresult, bootmean, bootbrd,
    bootrecover, jitterresult, jittermean, jitterbrd, jitterrecover,
    subsetresult, subsetmean, subsetbrd, subsetrecover, bojitresult,
    bojitmean, bojitbrd, bojitrecover, noiseresult, noisemean, 
    noisebrd, noiserecover}.

  \item{result}{clustering result; full output of the selected
    \code{clustermethod} for the original data set.}
  \item{partition}{partition parameter of the selected \code{clustermethod}
    (note that this is only meaningful for partitioning clustering methods).}
  \item{nc}{number of clusters in original data (including noise
    component if \code{noisemethod=TRUE}).}
  \item{nccl}{number of clusters in original data (not including noise
    component if \code{noisemethod=TRUE}).}
  \item{clustermethod, B, noisemethod, bootmethod, multipleboot, dissolution,
    recover}{input parameters, see above.}
  \item{bootresult}{matrix of Jaccard similarities for
    \code{bootmethod="boot"}. Rows correspond to clusters in the
    original data set. Columns correspond to bootstrap runs.}
  \item{bootmean}{clusterwise means of the \code{bootresult}.}
  \item{bootbrd}{clusterwise number of times a cluster has been dissolved.}
  \item{bootrecover}{clusterwise number of times a cluster has been
    successfully recovered.}
  \item{subsetresult, subsetmean, etc.}{same as \code{bootresult,
      bootmean, etc.}, but for the other resampling methods.} 
}
\references{
  Hennig, C. (2007) Cluster-wise assessment of cluster
  stability. \emph{Computational Statistics and Data Analysis},
  52, 258-271.

  Hennig, C. (2008)  Dissolution point and isolation robustness:
  robustness criteria for general cluster analysis methods.
  \emph{Journal of Multivariate Analysis} 99, 1154-1176.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\seealso{
  \code{\link{dist}},
  interface functions:
  \code{\link{kmeansCBI}}, \code{\link{hclustCBI}},
  \code{\link{hclusttreeCBI}}, \code{\link{disthclustCBI}},
  \code{\link{noisemclustCBI}}, \code{\link{distnoisemclustCBI}},
  \code{\link{claraCBI}}, \code{\link{pamkCBI}},
  \code{\link{trimkmeansCBI}}, \code{\link{disttrimkmeansCBI}},
  \code{\link{dbscanCBI}}, \code{\link{mahalCBI}}
}
\examples{  
  set.seed(20000)
  face <- rFace(50,dMoNo=2,dNoEy=0,p=2)
  cf1 <- clusterboot(face,B=3,bootmethod=
          c("boot","noise","jitter"),clustermethod=kmeansCBI,
          krange=5,seed=15555)
% For a serious application, choose a larger B!
  print(cf1)
  plot(cf1)
% cf1$result$result is the k-means clustering output for the original
% data; cf1$result is the output of bootkmeans on these data.
  cf2 <- clusterboot(dist(face),B=3,bootmethod=
          "subset",clustermethod=disthclustCBI,
          k=5, cut="number", method="average", showplots=TRUE, seed=15555)
  print(cf2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



