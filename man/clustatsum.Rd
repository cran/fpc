\name{clustatsum}
\alias{clustatsum}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compute and format cluster validation statistics}
\description{
  \code{clustatsum} computes cluster validation statistics by running
  \code{\link{cqcluster.stats}}, 
  and potentially \code{\link{distrsimilarity}}, and collecting some key
  statistics values with a somewhat different nomenclature.

  This was implemented as a helper function for use inside of
  \code{\link{clusterbenchstats}} and \code{\link{cgrestandard}}.
}
\usage{
clustatsum(datadist=NULL,clustering,noisecluster=FALSE,
                       datanp=NULL,npstats=FALSE,dnnk=2,
                       pamcrit=TRUE,...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{datadist}{distances on which validation-measures are based, \code{dist}
    object or distance matrix. If \code{NULL}, this is computed from
    \code{datanp}; at least one of \code{datadist} and \code{datanp}
    must be specified.}
  \item{clustering}{an integer vector of length of the number of cases,
    which indicates a clustering. The clusters have to be numbered
    from 1 to the number of clusters.}
  \item{noisecluster}{logical. If \code{TRUE}, it is assumed that the
    largest cluster number in \code{clustering} denotes a 'noise
    class', i.e. points that do not belong to any cluster. These points
    are not taken into account for the computation of all functions of
    within and between cluster distances including the validation
    indexes.} 
  \item{datanp}{optional observations times variables data matrix, see
    \code{npstats}.} 
  \item{npstats}{logical. If \code{TRUE}, \code{\link{distrsimilarity}}
    is called and the two statistics computed there are added to the
    output. These are based on \code{datanp} and require \code{datanp}
    to be specified.}
  \item{dnnk}{\code{nnk}-argument to be passed on to
    \code{\link{distrsimilarity}}.}
  \item{pamcrit}{\code{pamcrit}-argument to be passed on to
    \code{\link{cqcluster.stats}}.} 
  \item{...}{further arguments to be passed on to
    \code{\link{cqcluster.stats}}.}
}


\value{
  \code{clustatsum} returns a list. The components, as listed below, are
  outputs of \code{\link{summary.cquality}} with default parameters,
  which means that they are partly transformed versions of those given
  out by \code{\link{cqcluster.stats}}, i.e., their range is between 0
  and 1 and large values are good. Those from
  \code{\link{distrsimilarity}} are computed with
  \code{largeisgood=TRUE}, correspondingly.

  \item{avewithin}{average distance within clusters (reweighted so
  that every observation, rather than every distance, has the same weight).}
  \item{mnnd}{average distance to \code{nnk}th nearest neighbour within
    cluster.}
 \item{cvnnd}{coefficient of variation of dissimilarities to
   \code{nnk}th nearest wthin-cluster neighbour, measuring uniformity of
   within-cluster densities, weighted over all clusters, see Sec. 3.7 of
   Hennig (2017).}
  \item{maxdiameter}{maximum cluster diameter.}
  \item{widestgap}{widest within-cluster gap or average of cluster-wise
    widest within-cluster gap, depending on parameter \code{averagegap}.}
  \item{sindex}{separation index, see argument \code{sepindex}.}
  \item{minsep}{minimum cluster separation.}
  \item{asw}{average silhouette
    width. See \code{\link{silhouette}}.}
  \item{dindex}{this index measures to what extent the density decreases
    from the cluster mode to the outskirts; I-densdec in Sec. 3.6 of
    Hennig (2017); low values are good.}
  \item{denscut}{this index measures whether cluster boundaries run
    through density valleys; I-densbound in Sec. 3.6 of Hennig (2017); low
    values are good.}
  \item{highdgap}{this measures whether there is a large within-cluster
    gap with high density on both sides; I-highdgap in Sec. 3.6 of
    Hennig (2017); low values are good.}
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
  \item{pamc}{average distance to cluster centroid.}
  \item{kdnorm}{Kolmogorov distance between distribution of
    within-cluster Mahalanobis
    distances and appropriate chi-squared distribution, aggregated over
    clusters (I am grateful to Agustin Mayo-Iscar for the idea).}
  \item{kdunif}{Kolmogorov distance between distribution of distances to
    \code{nnk}th nearest within-cluster neighbor and appropriate
    Gamma-distribution, see Byers and Raftery (1998), aggregated over
    clusters.}
}
\references{

  Halkidi, M., Batistakis, Y., Vazirgiannis, M. (2001) On Clustering
  Validation Techniques, \emph{Journal of Intelligent Information
    Systems}, 17, 107-145.

  Hennig, C. (2017) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Proceedings of ASMDA 2017}, 501-520,
    \url{https://arxiv.org/abs/1703.09282}
  
    Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
  An Introduction to Cluster Analysis". Wiley, New York.
  
  Meila, M. (2007) Comparing clusterings?an information based distance,
  \emph{Journal of Multivariate Analysis}, 98, 873-895.
  
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{cqcluster.stats}}, \code{\link{distrsimilarity}}
}
\examples{  
  set.seed(20000)
  options(digits=3)
  face <- rFace(20,dMoNo=2,dNoEy=0,p=2)
  dface <- dist(face)
  complete3 <- cutree(hclust(dface),3)
  clustatsum(dface,complete3)
  
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}


