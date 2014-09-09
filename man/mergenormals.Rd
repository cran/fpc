\name{mergenormals}
\alias{mergenormals}
\alias{summary.mergenorm}
\alias{print.summary.mergenorm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Clustering by merging Gaussian mixture components}
\description{
  Clustering by merging Gaussian mixture components; computes all
  methods introduced in Hennig (2010) from an initial mclust
  clustering. See details section for details.
}
\usage{
  mergenormals(xdata, mclustsummary=NULL, 
                         clustering, probs, muarray, Sigmaarray, z,
                         method=NULL, cutoff=NULL, by=0.005,
                         numberstop=NULL, renumber=TRUE, M=50, ...)

  \method{summary}{mergenorm}(object, ...)

  \method{print}{summary.mergenorm}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{data (something that can be coerced into a matrix).}
  \item{mclustsummary}{output object from
    \code{\link[mclust]{summary.mclustBIC}} for \code{xdata}. Either
    \code{mclustsummary} or all of \code{clustering},
    \code{probs}, \code{muarray}, \code{Sigmaarray} and \code{z} need
    to be specified (the latter are obtained from \code{mclustsummary}
    if they are not provided). I am not aware of restrictions of the
    usage of
    \code{\link[mclust]{mclustBIC}} to produce an initial clustering;
    covariance matrix models can be restricted and a noise component can be
    included if desired, although I have probably not tested all
    possibilities.  
}
  \item{clustering}{vector of integers. Initial assignment of data to
    mixture components.}
  \item{probs}{vector of component proportions (for all components;
    should sum up to one).}
  \item{muarray}{matrix of component means (rows).}
  \item{Sigmaarray}{array of component covariance matrices (third
    dimension refers to component number).}
  \item{z}{matrix of observation- (row-)wise posterior probabilities of
    belonging to the components (columns).}
  \item{method}{one of \code{"bhat"}, \code{"ridge.uni"},
    \code{"ridge.ratio"}, \code{"demp"}, \code{"dipuni"},
    \code{"diptantrum"}, \code{"predictive"}. See details.}
  \item{cutoff}{numeric between 0 and 1. Tuning constant, see details
    and Hennig (2010). If not specified, the default values given in (9)
    in Hennig (2010) are used.}
  \item{by}{real between 0 and 1. Interval width for density computation
    along the ridgeline, used for methods \code{"ridge.uni"} and
    \code{"ridge.ratio"}. Methods \code{"dipuni"} and
    \code{"diptantrum"} require ridgeline computations and use it as well.} 
  \item{numberstop}{integer. If specified, \code{cutoff} is ignored and
    components are merged until the number of clusters specified here is
    reached.}
  \item{renumber}{logical. If \code{TRUE} merged clusters are renumbered
    from 1 to their number. If not, numbers of the original clustering
    are used (numbers of components that were merged into others then
    will not appear).}
  \item{M}{integer. Number of times the dataset is divided into two
    halves. Used if \code{method="predictive"}.} 
  \item{...}{additional optional parameters to pass on to
    \code{ridgeline.diagnosis} or \code{mixpredictive} (in
    \code{mergenormals}).}
  \item{object}{object of class \code{mergenorm}, output of
    \code{mergenormals}.}
  \item{x}{object of class \code{summary.mergenorm}, output of
    \code{summary.mergenorm}.} 
}

\value{
  \code{mergenormals} gives out an object of class \code{mergenorm},
  which is a List with components
  \item{clustering}{integer vector. Final clustering.}
  \item{clusternumbers}{vector of numbers of remaining clusters. These
    are given in terms of the original clusters even of
    \code{renumber=TRUE}, in which case they may be needed to understand
    the numbering of some further components, see below.}
  \item{defunct.components}{vector of numbers of components that were
    "merged away".}
  \item{valuemerged}{vector of values of the merging criterion (see
    details) at which components were merged.}
  \item{mergedtonumbers}{vector of numbers of clusters to which the
    original components were merged.}
  \item{parameters}{a list, if \code{mclustsummary} was provided. Entry
    no. i refers to number i in \code{clusternumbers}. The list entry i
    contains the parameters  of the original mixture components that
    make up cluster i, as extracted by
    \code{\link{extract.mixturepars}}.}
  \item{predvalues}{vector of prediction strength values for
    clusternumbers from 1 to the number of components in the original
    mixture, if \code{method=="predictive"}. See
    \code{\link{mixpredictive}}.}
  \item{orig.decisionmatrix}{square matrix with entries giving the
    original values of the merging criterion (see details) for every pair
    of original mixture components.}
  \item{new.decisionmatrix}{square matrix as \code{orig.decisionmatrix},
    but with final entries; numbering of rows and columns corresponds to
    \code{clusternumbers}; all entries corresponding to other rows and
    columns can be ignored.}
  \item{probs}{final cluster values of \code{probs} (see arguments)
    for merged components, generated by (potentially repeated) execution
    of  \code{\link{mergeparameters}} out of the original
    ones. Numbered according to \code{clusternumbers}.}
  \item{muarray}{final cluster means, analogous to \code{probs}.}
  \item{Sigmaarray}{final cluster covariance matrices, analogous to
    \code{probs}.}
  \item{z}{final matrix of posterior probabilities of observations
    belonging to the clusters, analogous to \code{probs}.}
  \item{noise}{logical. If \code{TRUE}, there was a noise component
    fitted in the initial mclust clustering (see help for
    \code{initialization} in \code{\link[mclust]{mclustBIC}}). In this
    case, a cluster number 0 indicates noise. noise is ignored by the
    merging methods and kept as it was originally.}
  \item{method}{as above.}  
  \item{cutoff}{as above.}
  
  \code{summary.mergenorm} gives out a list with components
  \code{clustering, clusternumbers, defunct.components, valuemerged,
    mergedtonumbers, predvalues, probs, muarray, Sigmaarray, z, noise,
    method, cutoff} as above, plus \code{onc} (original number of
    components) and \code{mnc} (number of clusters after merging).

}

\details{
  Mixture components are merged in a hierarchical fashion. The merging
  criterion is computed for all pairs of current clusters and the two
  clusters with the highest criterion value (lowest, respectively, for
  \code{method="predictive"}) are merged. Then criterion values are
  recomputed for the merged cluster. Merging is continued until the
  criterion value to merge is below (or above, for
  \code{method="predictive"}) the cutoff value. Details are given in
  Hennig (2010). The following criteria are offered, specified by the
  \code{method}-argument.
  \describe{
    \item{"ridge.uni"}{components are only merged if their mixture is
      unimodal according to Ray and Lindsay's (2005) ridgeline theory,
      see \code{\link{ridgeline.diagnosis}}. This ignores argument
      \code{cutoff}.}
    \item{"ridge.ratio"}{ratio between density minimum between
      components and minimum of density maxima according to Ray and
      Lindsay's (2005) ridgeline theory, see
      \code{\link{ridgeline.diagnosis}}. }      
    \item{"bhat"}{Bhattacharyya upper bound on misclassification
      probability between two components, see
      \code{\link{bhattacharyya.matrix}}.}
    \item{"demp"}{direct estimation of misclassification probability
      between components, see Hennig (2010).}
    \item{"dipuni"}{this uses \code{method="ridge.ratio"} to decide
      which clusters to merge but stops merging according to the p-value of
      the dip test computed as in Hartigan and Hartigan (1985), see
      \code{\link[diptest]{dip.test}}.}
    \item{"diptantrum"}{as \code{"dipuni"}, but p-value of dip test
      computed as in Tantrum, Murua and Stuetzle (2003), see
      \code{\link{dipp.tantrum}}.}
    \item{"predictive"}{this uses \code{method="demp"} to decide which
      clusters to merge but stops merging according to the value of
      prediction strength (Tibshirani and Walther, 2005) as computed in
      \code{\link{mixpredictive}}.}
    }
}

  \references{
J. A. Hartigan and P. M. Hartigan (1985) The Dip Test of
     Unimodality, \emph{Annals of Statistics}, 13, 70-84.

Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
  
Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.


Tantrum, J., Murua, A. and Stuetzle, W. (2003) Assessment and 
Pruning of Hierarchical Model Based Clustering, \emph{Proceedings of the 
ninth ACM SIGKDD international conference on Knowledge discovery and 
data mining}, Washington, D.C., 197-205.

Tibshirani, R. and Walther, G. (2005) Cluster Validation by 
Prediction Strength, \emph{Journal of Computational and Graphical 
  Statistics}, 14, 511-528.

}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\examples{
  require(mclust)
  require(MASS)
  options(digits=3)
  data(crabs)
  dc <- crabs[,4:8]
  cm <- mclustBIC(crabs[,4:8],G=9,modelNames="EEE")
  scm <- summary(cm,crabs[,4:8])
  cmnbhat <- mergenormals(crabs[,4:8],scm,method="bhat")
  summary(cmnbhat)
  cmndemp <- mergenormals(crabs[,4:8],scm,method="demp")
  summary(cmndemp)
# Other methods take a bit longer, but try them!
# The values of by and M below are still chosen for reasonably fast execution.
# cmnrr <- mergenormals(crabs[,4:8],scm,method="ridge.ratio",by=0.05)
# cmd <- mergenormals(crabs[,4:8],scm,method="dip.tantrum",by=0.05)
# cmp <- mergenormals(crabs[,4:8],scm,method="predictive",M=3)
}
\keyword{multivariate}
\keyword{cluster}


