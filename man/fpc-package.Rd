\name{fpc-package}
\alias{fpc-package}
%- Also NEED an `\alias' for EACH other topic documented here.
\docType{package}
\title{fpc package overview}
\description{
  Here is a list of the main functions in package fpc. Most other
  functions are auxiliary functions for these.
}
\section{Clustering methods}{
  \describe{
    \item{dbscan}{Computes DBSCAN density based clustering as introduced
      in Ester et al. (1996).}
    \item{fixmahal}{Mahalanobis Fixed Point Clustering, Hennig and
      Christlieb (2002), Hennig (2005).}
    \item{fixreg}{Regression Fixed Point Clustering, Hennig (2003).}
    \item{flexmixedruns}{This fits a latent class model to
      data with mixed type continuous/nominal variables. Actually it
      calls a method for \code{\link[flexmix]{flexmix}}.}
    \item{mergenormals}{Clustering by merging components of a Gaussian
      mixture, see Hennig (2010).}
    \item{regmix}{ML-fit of a mixture of linear regression models, see
      DeSarbo and Cron (1988).}
    }}
\section{Cluster validity indexes and estimation of the number of clusters}{
  \describe{
    \item{cluster.stats}{This computes several cluster validity
      statistics from a clustering and a dissimilarity matrix including
      the Calinski-Harabasz index, the adjusted Rand index and other
      statistics explained in Gordon (1999) as well as several
      characterising
      measures such as average between cluster and within cluster
      dissimilarity and separation. See also \code{\link{calinhara}},
      \code{\link{dudahart2}} for specific indexes, and a new version
      \code{\link{cqcluster.stats}} that computes some more indexes and
      statistics used for computing them. There's also
      \code{\link{distrsimilarity}}, which computes within-cluster
      dissimilarity to the Gaussian and uniform distribution.}
    \item{prediction.strength}{Estimates the number of clusters by
      computing the prediction strength of a
      clustering of a dataset into different numbers of components for
      various clustering methods, see
      Tibshirani and Walther (2005). In fact, this is more flexible than
      what is in the original paper, because it can use
      point classification schemes that work better with clustering
      methods other than k-means.}
    \item{nselectboot}{Estimates the number of clusters by bootstrap
      stability selection, see Fang and Wang (2012). This is quite
      flexible regarding clustering methods and point classification
      schemes and also allows for dissimilarity data.}
    \item{clusterbenchstats}{This runs many clustering methods (to be
      specifed by the user) with many numbers of clusters on a dataset
      and produces standardised and comparable versions of many cluster
      validity indexes (see Hennig 2019, Akhanli and Hennig 2020).
      This is done by means of
      producing random clusterings on the given data, see
      \code{\link{stupidkcentroids}} and \code{\link{stupidknn}}. It
      allows to compare many 
      clusterings based on many different potential desirable features
      of a clustering. \code{\link{print.valstat}} allows to compute an
      aggregated index with user-specified weights.}
    }}
\section{Cluster visualisation and validation}{
  \describe{
    \item{clucols}{Sets of colours and symbols useful for cluster plotting.}
    \item{clusterboot}{Cluster-wise stability assessment of a
      clustering. Clusterings are performed on resampled data to see for
      every cluster of the original dataset how well this is
      reproduced. See Hennig (2007) for details.}
    \item{cluster.varstats}{Extracts variable-wise information for every
      cluster in order to help with cluster interpretation.}
    \item{plotcluster}{Visualisation of a clustering or grouping in data
      by various linear projection methods that optimise the separation
      between clusters, or between a single cluster and the rest of the
      data according to Hennig (2004) including classical methods such
      as discriminant coordinates. This calls the function
      \code{\link{discrproj}}, which is a bit more flexible but doesn't
      produce a plot itself.}
    \item{ridgeline.diagnosis}{Plots and diagnostics for assessing
      modality of Gaussian mixtures, see Ray and Lindsay (2005).}
    \item{weightplots}{Plots to diagnose component separation in
      Gaussian mixtures, see Hennig (2010).}
    \item{localshape}{Local shape matrix, can be used for finding
      clusters in connection with function \code{ics} in package
  \code{ICS}, see Hennig's
  discussion and rejoinder of Tyler et al. (2009).}
    }}
\section{Useful wrapper functions for clustering methods}{
  \describe{
    \item{kmeansCBI}{This and other "CBI"-functions (see the
      \code{\link{kmeansCBI}}-help page) are unified wrappers for
      various clustering methods in R that may be useful because they do
      in one step for what you normally may need to do a bit more in R
      (for example fitting a Gaussian mixture with noise component in
      package mclust).}
    \item{kmeansruns}{This calls \code{\link{kmeans}} for the k-means
      clustering method and includes estimation of the number of
      clusters and finding an optimal solution from several starting
      points.}
    \item{pamk}{This calls \code{\link[cluster]{pam}} and
      \code{\link[cluster]{clara}} for the partitioning around medoids 
      clustering method (Kaufman and Rouseeuw, 1990) and includes two
      different ways of estimating the number of clusters.}
    }}

\references{
Akhanli, S. and Hennig, C. (2020) Calibrating and aggregating cluster
validity indexes for context-adapted comparison of clusterings.
\emph{Statistics and Computing}, 30, 1523-1544,
\url{https://link.springer.com/article/10.1007/s11222-020-09958-2},
\url{https://arxiv.org/abs/2002.01822}

  DeSarbo, W. S. and Cron, W. L. (1988) A maximum likelihood methodology
  for clusterwise linear regression, \emph{Journal of
    Classification} 5, 249-282.
  
  Ester, M., Kriegel, H.-P., Sander, J. and Xu, X.
(1996). A Density-Based Algorithm for Discovering Clusters in Large Spatial
Databases with Noise.
\emph{Proceedings of 2nd International Conference on
  Knowledge Discovery and Data
Mining (KDD-96).}

  Fang, Y. and Wang, J. (2012) Selection of the number of clusters via
  the bootstrap method. \emph{Computational Statistics and Data
    Analysis}, 56, 468-477.


Gordon, A. D. (1999) \emph{Classification}, 2nd ed. Chapman and Hall.

Hennig, C. (2003) Clusters, outliers and regression:
  fixed point clusters, \emph{Journal of
    Multivariate Analysis} 86, 183-212.
  
  Hennig, C. (2004) Asymmetric linear dimension reduction for classification.
  \emph{Journal of Computational and Graphical Statistics}, 13, 930-945 .

   Hennig, C. (2005) Fuzzy and Crisp Mahalanobis Fixed Point Clusters,
  in Baier, D., Decker, R., and Schmidt-Thieme, L. (eds.):
  \emph{Data Analysis and Decision Support}. Springer, Heidelberg,
  47-56.
  
    Hennig, C. (2007) Cluster-wise assessment of cluster
  stability. \emph{Computational Statistics and Data Analysis},
  52, 258-271.

Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.

Hennig, C. (2019) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Data Analysis and Applications 1: Clustering and Regression,
  Modeling-estimating, Forecasting and Data Mining, Volume 2}, Wiley,
  New York 1-24,
    \url{https://arxiv.org/abs/1703.09282}

  
  Hennig, C. and Christlieb, N. (2002) Validating visual clusters in
  large datasets: Fixed point clusters of spectral features,
  \emph{Computational Statistics and Data Analysis} 40, 723-739.

  Kaufman, L. and Rousseeuw, P.J. (1990). "Finding Groups in Data:
  An Introduction to Cluster Analysis". Wiley, New York.

  Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.

  Tibshirani, R. and Walther, G. (2005) Cluster Validation by 
Prediction Strength, \emph{Journal of Computational and Graphical 
  Statistics}, 14, 511-528.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}



