\name{fixreg}
\alias{fixreg}
\alias{summary.rfpc}
\alias{plot.rfpc}
\alias{fpclusters.rfpc}
\alias{print.summary.rfpc}
\alias{print.rfpc}
\alias{rfpi}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Linear Regression Fixed Point Clusters}
\description{
  Computes linear regression fixed point clusters (FPCs), i.e.,
  subsets of the data, which consist exactly of the
  non-outliers w.r.t. themselves, and may be interpreted
  as generated from a homogeneous linear regression relation
  between independent and dependent variable.  FPCs may
  overlap, are not necessarily exhausting and
  do not need a specification of the number of clusters.

  Note that while \code{fixreg} has lots of parameters, only one (or
  few) of them have usually to be specified, cf. the examples. The
  philosophy is to allow much flexibility, but to always provide 
  sensible defaults.
  
}
\usage{
fixreg(indep=rep(1,n), dep, n=length(dep),
                    p=ncol(as.matrix(indep)),
                    ca=NA, mnc=NA, mtf=3, ir=NA, irnc=NA,
                    irprob=0.95, mncprob=0.5, maxir=20000, maxit=5*n,
                    distcut=0.85, init.group=list(), 
                    ind.storage=FALSE, countmode=100, 
                    plot=FALSE)

\method{summary}{rfpc}(object, ...)

\method{print}{summary.rfpc}(x, maxnc=30, ...)

\method{plot}{rfpc}(x, indep=rep(1,n), dep, no, bw=TRUE,
                      main=c("Representative FPC No. ",no),
                      xlab="Linear combination of independents",
                      ylab=deparse(substitute(indep)),
                      xlim=NULL, ylim=range(dep), 
                      pch=NULL, col=NULL,...)

\method{fpclusters}{rfpc}(object, indep=NA, dep=NA, ca=object$ca, ...)

rfpi(indep, dep, p, gv, ca, maxit, plot) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{indep}{numerical matrix or vector. Independent
    variables.
    Leave out for clustering one-dimensional data.
    \code{fpclusters.rfpc}
    does not need specification of \code{indep} if \code{fixreg}
    was run with \code{ind.storage=TRUE}.}
  \item{dep}{numerical vector. Dependent variable.
    \code{fpclusters.rfpc}
    does not need specification of \code{dep} if \code{fixreg}
    was run with \code{ind.storage=TRUE}.}
  \item{n}{optional positive integer.
    Number of cases.}
  \item{p}{optional positive integer.
    Number of independent variables.}
  \item{ca}{optional positive number. Tuning constant, specifying
    required cluster
    separation. By default determined automatically as a
    function of \code{n} and \code{p}, see function \code{\link{can}},
    Hennig (2002a).}
  \item{mnc}{optional positive integer. Minimum size of clusters
    to be reported.
    By default determined automatically as a function of
    \code{mncprob}. See Hennig (2002a).}
  \item{mtf}{optional positive integer. FPCs must be found at
    least \code{mtf} times to be reported by \code{summary.rfpc}.}
  \item{ir}{optional positive integer. Number of algorithm runs.
    By default determined
    automatically as a function of \code{n}, \code{p}, \code{irnc},
    \code{irprob}, \code{mtf},
    \code{maxir}. See function \code{\link{itnumber}} and Hennig (2002a).}
  \item{irnc}{optional positive integer. Size of the smallest
    cluster to be found with
    approximated probability \code{irprob}.}
  \item{irprob}{optional value between 0 and 1. Approximated
    probability for a cluster of size \code{irnc} to be found.}
  \item{mncprob}{optional value between 0 amd 1. Approximated
    probability for a cluster of size \code{mnc} to be found.}
  \item{maxir}{optional integer. Maximum number of algorithm runs.}
  \item{maxit}{optional integer. Maximum number of iterations
    per algorithm run (usually an FPC is found much earlier).}
  \item{distcut}{optional value between 0 and 1. A similarity
    measure between FPCs, given in Hennig (2002a), and the corresponding
    Single Linkage groups of FPCs with similarity larger
    than \code{distcut} are computed.
    A single representative FPC is selected for each group.}
  \item{init.group}{optional list of logical vectors of length
    \code{n}.
    Every vector indicates a starting configuration for the fixed
    point algorithm. This can be used for datasets with high
    dimension, where the vectors of \code{init.group} indicate cluster
    candidates found by graphical inspection or background
    knowledge.}
  \item{ind.storage}{optional logical. If \code{TRUE},
    then all indicator
    vectors of found FPCs are given in the value of \code{fixreg}.
    May need lots of memory, but is a bit faster.}
  \item{countmode}{optional positive integer. Every \code{countmode}
    algorithm runs \code{fixreg} shows a message.}
  \item{plot}{optional logical. If \code{TRUE}, you get a scatterplot
    of first independent vs. dependent variable at each iteration.}
  \item{object}{object of class \code{rfpc}, output of \code{fixreg}.}
  \item{x}{object of class \code{rfpc}, output of \code{fixreg}.}
  \item{maxnc}{positive integer. Maximum number of FPCs
    to be reported.}
  \item{no}{positive integer. Number of the representative FPC to
    be plotted.}
  \item{bw}{optional logical. If \code{TRUE}, plot is black/white,
    FPC is
    indicated by different symbol. Else FPC is indicated red.}
  \item{main}{plot title.}
  \item{xlab}{label for x-axis.}
  \item{ylab}{label for y-axis.}
  \item{xlim}{plotted range of x-axis. If \code{NULL}, the range of the
    plotted linear combination of independent variables is used.}
  \item{ylim}{plotted range of y-axis.}
  \item{pch}{plotting symbol, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{col}{plotting color, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{gv}{logical vector of length \code{n}. Indicates the initial
    configuration for the fixed point algorithm.} 
  \item{...}{additional parameters to be passed to \code{plot}
    (no effects elsewhere).}
}
\details{A linear regression FPC is a data subset 
  that reproduces itself under the following operation: \cr
  Compute linear regression and error variance estimator for the data
  subset, and compute all points of the dataset for which the squared
  residual is smaller than \code{ca} times the error variance.\cr
  Fixed points of this operation can be considered as clusters,
  because they contain only
  non-outliers (as defined by the above mentioned procedure) and all other
  points are outliers w.r.t. the subset. \cr
  \code{fixreg} performs \code{ir} fixed point algorithms started from
  random subsets of size \code{p+2} to look for
  FPCs. Additionally an algorithm is started from the whole dataset,
  and algorithms are started from the subsets specified in
  \code{init.group}. \cr
  Usually some of the FPCs are unstable, and more than one FPC may
  correspond to the same significant pattern in the data. Therefore the
  number of FPCs is reduced: FPCs with less than \code{mnc} points are
  ignored. Then a similarity matrix is computed between the remaining
  FPCs. Similarity between sets is defined as the ratio between
  2 times size of
  intersection and the sum of sizes of both sets. The Single Linkage
  clusters (groups)
  of level \code{distcut} are computed, i.e. the connectivity
  components of the graph where edges are drawn between FPCs with
  similarity larger than \code{distcut}. Groups of FPCs whose members
  are found \code{mtf} times or more are considered as stable enough.
  A representative FPC is
  chosen for every Single Linkage cluster of FPCs according to the
  maximum expectation ratio \code{ser}. \code{ser} is the ratio between
  the number of findings of an FPC and the estimated
  expectation of the number of findings of an FPC of this size,
  called \emph{expectation ratio} and
  computed by \code{\link{clusexpect}}.\cr
  Usually only the representative FPCs of stable groups
  are of interest. \cr
  The choice of the involved tuning constants such as \code{ca} and
  \code{ir} is discussed in detail in Hennig (2002a). Statistical theory
  is presented in Hennig (2003).\cr
  Generally, the default settings are recommended for
  \code{fixreg}. In cases where they lead to a too large number of
  algorithm runs (e.g., always for \code{p>4}), the use of
  \code{init.group} together with \code{mtf=1} and \code{ir=0}
  is useful. Occasionally, \code{irnc} may be chosen
  smaller than the default,
  if smaller clusters are of interest, but this may lead to too many
  clusters and too many algorithm runs. Decrease of
  \code{ca} will often lead to too many clusters, even for homogeneous
  data. Increase of \code{ca} will produce only very strongly
  separated clusters. Both may be of interest occasionally.

  \code{rfpi} is called by \code{fixreg} for a single fixed point
  algorithm and will usually not be executed alone.

  \code{summary.rfpc} gives a summary about the representative FPCs of
  stable groups.

  \code{plot.rfpc} is a plot method for the representative FPC of stable
  group 
  no. \code{no}. It produces a scatterplot of the linear combination of
  independent variables determined by the regression coefficients of the
  FPC vs. the dependent variable. The regression line and the region
  of non-outliers determined by \code{ca} are plotted as well.

  \code{fpclusters.rfpc} produces a list of indicator vectors for the
  representative FPCs of stable groups.
}
\value{
  \code{fixreg} returns an object of class \code{rfpc}. This is a list
  containing the components \code{nc, g, coefs, vars, nfound, er, tsc,
    ncoll, grto, imatrix, smatrix, stn, stfound, sfpc, ssig, sto, struc,
    n, p, ca, ir, mnc, mtf, distcut}.

  \code{summary.rfpc} returns an object of class \code{summary.rfpc}.
  This is a list containing the components \code{coefs, vars, stfound,
    stn, sn, ser, tsc, sim, ca, ir, mnc, mtf}.

  \code{fpclusters.rfpc} returns a list of indicator vectors for the
  representative FPCs of stable groups.

  \code{rfpi} returns a list with the components \code{coef, var, g,
    coll, ca}.
  
  \item{nc}{integer. Number of FPCs.}
  \item{g}{list of logical vectors. Indicator vectors of FPCs. \code{FALSE}
    if \code{ind.storage=FALSE}.}
  \item{coefs}{list of numerical vectors. Regression coefficients of
    FPCs. In \code{summary.rfpc}, only for representative
    FPCs of stable groups and sorted according to
    \code{stfound}.}
  \item{vars}{list of numbers. Error variances of FPCs. In
    \code{summary.rfpc}, only for representative
    FPCs of stable groups and sorted according to
    \code{stfound}.}
  \item{nfound}{vector of integers. Number of findings for the FPCs.}
  \item{er}{numerical vector. Expectation ratios of FPCs. Can be
    taken as a stability measure.}
  \item{tsc}{integer. Number of algorithm runs leading to too small or
    too seldom found FPCs.}
  \item{ncoll}{integer. Number of algorithm runs where collinear
    regressor matrices occurred.}
  \item{grto}{vector of integers. Numbers of FPCs to which algorithm
    runs led, which were started by \code{init.group}.}
  \item{imatrix}{vector of integers. Size of intersection between
    FPCs. See \code{\link{sseg}}.}
  \item{smatrix}{numerical vector. Similarities between
    FPCs. See \code{\link{sseg}}.}
  \item{stn}{integer. Number of representative FPCs of stable groups. In
    \code{summary.rfpc} sorted according to \code{stfound}.}
  \item{stfound}{vector of integers. Number of findings of members of
    all groups of FPCs. In
    \code{summary.rfpc} sorted according to \code{stfound}.}
  \item{sfpc}{vector of integers. Numbers of representative FPCs.}
  \item{ssig}{vector of integers. As \code{sfpc}, but only for stable
    groups.}
  \item{sto}{vector of integers. Number of representative FPC of most,
    2nd most, ..., often found group of FPCs.}
  \item{struc}{vector of integers. Number of group an FPC belongs to.}
  \item{n}{see arguments.}
  \item{p}{see arguments.}
  \item{ca}{see arguments.}
  \item{ir}{see arguments.}
  \item{mnc}{see arguments.}
  \item{mtf}{see arguments.}
  \item{distcut}{see arguments.}

  \item{sn}{vector of integers. Number of points of representative
    FPCs.}
  \item{ser}{numerical vector. Expectation ratio for stable groups.}
  \item{sim}{vector of integers. Size of intersections between
    representative FPCs of stable groups. See \code{\link{sseg}}.}

  \item{coef}{vector of regression coefficients.}
  \item{var}{error variance.}
  \item{g}{logical indicator vector of iterated FPC.}
  \item{coll}{logical. \code{TRUE} means that singular covariance
    matrices occurred during the iterations.}
}
\references{
  Hennig, C. (2002) Fixed point clusters for linear regression:
  computation and comparison, \emph{Journal of
    Classification} 19, 249-276.

  Hennig, C. (2003) Clusters, outliers and regression:
  fixed point clusters, \emph{Journal of
    Multivariate Analysis} 86, 183-212.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}
\seealso{
\code{\link{fixmahal}} for fixed point clusters in the usual setup
(non-regression).

\code{\link{regmix}} for clusterwiese linear regression by mixture
modeling ML.

\code{\link{can}}, \code{\link{itnumber}} for computation of the default
settings.  

\code{\link{clusexpect}} for estimation of the expected number of
findings of an FPC of given size.

\code{\link{itnumber}} for the generation of the number of fixed point
algorithms.

\code{\link{minsize}} for the smallest FPC size to be found with a given
probability..

\code{\link{sseg}} for indexing the similarity/intersection vectors
computed by \code{fixreg}.
}
\examples{
set.seed(190000)
data(tonedata)
# Note: If you do not use the installed package, replace this by
# tonedata <- read.table("(path/)tonedata.txt", header=TRUE)
attach(tonedata)
tonefix <- fixreg(stretchratio,tuned,mtf=1,ir=20)
summary(tonefix)
# This is designed to have a fast example; default setting would be better.
# If you want to see more (and you have a bit more time),
# try out the following:
# set.seed(1000)
# tonefix <- fixreg(stretchratio,tuned)
## Default - good for these data
# summary(tonefix)
# plot(tonefix,stretchratio,tuned,1)
# plot(tonefix,stretchratio,tuned,2)
# plot(tonefix,stretchratio,tuned,3,bw=FALSE,pch=5) 
# toneclus <- fpclusters(tonefix,stretchratio,tuned)
# plot(stretchratio,tuned,col=1+toneclus[[2]])
# tonefix2 <- fixreg(stretchratio,tuned,distcut=1,mtf=1,countmode=50)
## Every found fixed point cluster is reported,
## no matter how instable it may be.
# summary(tonefix2)
# tonefix3 <- fixreg(stretchratio,tuned,ca=7)
## ca defaults to 10.07 for these data.
# summary(tonefix3)
# subset <- c(rep(FALSE,5),rep(TRUE,24),rep(FALSE,121))
# tonefix4 <- fixreg(stretchratio,tuned,
#                    mtf=1,ir=0,init.group=list(subset))
# summary(tonefix4)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{robust}% __ONLY ONE__ keyword per line
\keyword{regression}

