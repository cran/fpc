\name{valstat.object}
\alias{valstat.object}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cluster validation statistics - object}
\description{
  The objects of class \code{"valstat"} store cluster validation
  statistics from various clustering methods run with various numbers of
  clusters.
}
\section{GENERATION}{
  These objects are generated as part of the
  \code{\link{clusterbenchstats}}-output.
}

\section{METHODS}{
  The \code{valstat} class has methods for the following generic functions:
  \code{print}, \code{plot}, see \code{\link{plot.valstat}}.
}
\value{
  A legitimate \code{valstat} object is a list. The format of the list
  relies on the number of involved clustering methods, \code{nmethods},
  say, i.e., the length
  of the \code{method}-component explained below. The first
  \code{nmethods} elements of the \code{valstat}-list are just
  numbered. These are themselves lists that are numbered between 1 and
  the \code{maxG}-component defined below. Element \code{[[i]][[j]]}
  refers to the clustering from clustering method number i with number
  of clusters j. Every such element is a list
  with components 
  \code{avewithin, mnnd, cvnnd, maxdiameter, widestgap, sindex, minsep,
  asw, dindex, denscut, highdgap, pearsongamma, withinss, entropy}:
  Further optional components are \code{pamc, kdnorm, kdunif,
    dmode, aggregated}. All these are cluster validation indexes, as
  follows.

  \item{avewithin}{average distance within clusters (reweighted so
  that every observation, rather than every distance, has the same weight).}
  \item{mnnd}{average distance to \code{nnk}th nearest neighbour within
    cluster. (\code{nnk} is a parameter of
    \code{\link{cqcluster.stats}}, default 2.)}
  \item{cvnnd}{coefficient of variation of dissimilarities to
    \code{nnk}th nearest wthin-cluster neighbour, measuring uniformity of
    within-cluster densities, weighted over all clusters, see Sec. 3.7 of
    Hennig (2019). (\code{nnk} is a parameter of
    \code{\link{cqcluster.stats}}, default 2.)}
  \item{maxdiameter}{maximum cluster diameter.}
  \item{widestgap}{widest within-cluster gap or average of cluster-wise
    widest within-cluster gap, depending on parameter \code{averagegap}
    of \code{\link{cqcluster.stats}}, default \code{FALSE}.}
 \item{sindex}{separation index. Defined based on the distances for
    every point to the
    closest point not in the same cluster. The separation index is then
    the mean of the smallest proportion \code{sepprob} (parameter of
    \code{\link{cqcluster.stats}}, default 0.1) of these. See Hennig (2019).}
  \item{minsep}{minimum cluster separation.}
  \item{asw}{average silhouette
    width. See \code{\link{silhouette}}.}
  \item{dindex}{this index measures to what extent the density decreases
    from the cluster mode to the outskirts; I-densdec in Sec. 3.6 of
    Hennig (2019); low values are good.}
  \item{denscut}{this index measures whether cluster boundaries run
    through density valleys; I-densbound in Sec. 3.6 of Hennig (2019); low
    values are good.}
  \item{highdgap}{this measures whether there is a large within-cluster
    gap with high density on both sides; I-highdgap in Sec. 3.6 of
    Hennig (2019); low values are good.}
  \item{pearsongamma}{correlation between distances and a
    0-1-vector where 0 means same cluster, 1 means different clusters.
    "Normalized gamma" in Halkidi et al. (2001).}
  \item{withinss}{a generalisation of the within clusters sum
    of squares (k-means objective function), which is obtained if
    \code{d} is a Euclidean distance matrix.  For general distance
    measures, this is half
    the sum of the within cluster squared dissimilarities divided by the
    cluster size.}
  \item{entropy}{entropy of the distribution of cluster memberships,
    see Meila(2007).}
  \item{pamc}{average distance to cluster centroid, which is the
    observation that minimises this average distance.}
  \item{kdnorm}{Kolmogorov distance between distribution of
    within-cluster Mahalanobis
    distances and appropriate chi-squared distribution, aggregated over
    clusters (I am grateful to Agustin Mayo-Iscar for the idea).}
  \item{kdunif}{Kolmogorov distance between distribution of distances to
    \code{dnnk}th nearest within-cluster neighbor and appropriate
    Gamma-distribution, see Byers and Raftery (1998), aggregated over
    clusters. \code{dnnk} is parameter \code{nnk} of
    \code{\link{distrsimilarity}}, corresponding to \code{dnnk} of
    \code{\link{clusterbenchstats}}.}
  \item{dmode}{aggregated density mode index equal to
    \code{0.75*dindex+0.25*highdgap} before standardisation.}

  Furthermore, a \code{valstat} object
  has the following list components:   

  \item{maxG}{maximum number of clusters.}
  \item{minG}{minimum number of clusters (list entries below that number
    are empty lists).}
  \item{method}{vector of names (character strings) of clustering
    CBI-functions, see \code{\link{kmeansCBI}}.}
  \item{name}{vector of names (character strings) of clustering
    methods. These can be user-chosen names (see argument
    \code{methodsnames} in \code{\link{clusterbenchstats}}) and may
    distinguish different methods run by the same CBI-function but with
    different parameter values such as complete and average linkage for
    \code{\link{hclustCBI}}.}
  \item{statistics}{vector of names (character strings) of cluster
    validation indexes.} 

  
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
  \code{\link{clusterbenchstats}},
  \code{\link{plot.valstat}}.
}

\keyword{cluster}% at least one, from doc/KEYWORDS



