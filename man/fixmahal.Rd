\name{fixmahal}
\alias{fixmahal}
\alias{summary.mfpc}
\alias{plot.mfpc}
\alias{fpclusters.mfpc}
\alias{print.summary.mfpc}
\alias{print.mfpc}
\alias{fpmi}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mahalanobis Fixed Point Clusters}
\description{
  Computes Mahalanobis fixed point clusters (FPCs), i.e.,
  subsets of the data, which consist exactly of the
  non-outliers w.r.t. themselves, and may be interpreted
  as generated from a homogeneous normal population.
  FPCs may
  overlap, are not necessarily exhausting and
  do not need a specification of the number of clusters.

  Note that while \code{fixmahal} has lots of parameters, only one (or
  few) of them have usually to be specified, cf. the examples. The
  philosophy is to allow much flexibility, but to always provide 
  sensible defaults.
}
\usage{
fixmahal(dat, n = nrow(as.matrix(dat)), p = ncol(as.matrix(dat)), 
                      method = "fuzzy", cgen = "fixed",
                      ca = NA, ca2 = NA,
                      calpha = ifelse(method=="fuzzy",0.95,0.99),
                      calpha2 = 0.995,
                      pointit = TRUE, subset = n,
                      nc1 = 100+20*p,
                      startn = 18+p, mnc = floor(startn/2), 
                      mer = ifelse(pointit,0.1,0), 
                      distcut = 0.85, maxit = 5*n, iter = n*1e-5,
                      init.group = list(), 
                      ind.storage = TRUE, countmode = 100, 
                      plot = "none")


\method{summary}{mfpc}(object, ...)

\method{print}{summary.mfpc}(x, maxnc=30, ...)

\method{plot}{mfpc}(x, dat, no, bw=FALSE, main=c("Representative FPC No. ",no),
                    xlab=NULL, ylab=NULL,
                    pch=NULL, col=NULL, ...)

\method{fpclusters}{mfpc}(object, dat=NA, ca=object$ca, p=object$p, ...)

fpmi(dat, n = nrow(as.matrix(dat)), p = ncol(as.matrix(dat)),
                  gv, ca, ca2, method = "ml", plot,
                  maxit = 5*n, iter = n*1e-6) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{something that can be coerced to a
    numerical matrix or vector. Data matrix, rows are points, columns
    are variables.
    \code{fpclusters.rfpc}
    does not need specification of \code{dat} if \code{fixmahal}
    has been run with \code{ind.storage=TRUE}.}
  \item{n}{optional positive integer.
    Number of cases.}
  \item{p}{optional positive integer.
    Number of independent variables.}
  \item{method}{a string. \code{method="classical"} means 0-1 weighting
    of observations by Mahalanobis distances and use of the classical
    normal covariance estimator. \code{method="ml"} uses the
    ML-covariance estimator (division by \code{n} instead of \code{n-1}) 
    This is used in Hennig and Christlieb (2002).
    \code{method} can also be \code{"mcd"} or \code{"mve"},
    to enforce the use of robust centers and covariance matrices, see
    \code{\link{cov.rob}}. This is experimental, not recommended at the
    moment, may be very slowly and requires library \code{lqs}.
    The default is
    \code{method="fuzzy"}, where weighted means and covariance matrices
    are used (Hennig, 2005).
    The weights are computed by \code{\link{wfu}}, i.e., a
    function that is constant 1 for arguments smaller than \code{ca}, 0 for
    arguments larger than \code{ca2} and continuously linear in between.
    Convergence is only proven for \code{method="ml"} up to now.}  
  \item{cgen}{optional string. \code{"fixed"} means that the same tuning
    constant \code{ca} is used for all iterations. \code{"auto"} means
    that \code{ca} is generated dependently on the size of the current data
    subset in each iteration by \code{\link{cmahal}}. This is 
    experimental.}
  \item{ca}{optional positive number. Tuning constant, specifying
    required cluster
    separation. By default determined as \code{calpha}-quantile of the
    chisquared distribution with \code{p} degrees of freedom.}
  \item{ca2}{optional positive number. Second tuning constant needed if
    \code{method="fuzzy"}.
    By default determined as \code{calpha2}-quantile of the
    chisquared distribution with \code{p} degrees of freedom.}
  \item{calpha}{number between 0 and 1. See \code{ca}.}
  \item{calpha2}{number between 0 and 1, larger than \code{calpha}.
    See \code{ca2}.}
  \item{pointit}{optional logical. If \code{TRUE}, \code{subset} fixed point
    algorithms are started from initial configurations, which are built
    around single points of the dataset, cf. \code{\link{mahalconf}}.
    Otherwise, initial configurations are only specified by
    \code{init.group}.}
  \item{subset}{optional positive integer smaller or equal than \code{n}.
    Initial configurations for the fixed point algorithm
    (cf. \code{\link{mahalconf}}) are built from
    a subset of \code{subset} points from the data. No effect if
    \code{pointit=FALSE}. Default: all points.}
  \item{nc1}{optional positive integer. Tuning constant needed by
    \code{\link{cmahal}} to generate \code{ca} automatically. Only 
    needed for \code{cgen="auto"}.} 
  \item{startn}{optional positive integer. Size of the initial
    configurations. The default value is chosen to prevent that small 
    meaningless FPCs are found, but it should be decreased if 
    clusters of size smaller than the default value are of interest.} 
  \item{mnc}{optional positive integer. Minimum size of clusters
    to be reported.}
  \item{mer}{optional nonnegative number. FPCs (groups of them,
    respectively, see details)
    are only reported as stable if the ratio 
    of the number of their
    findings to their number of points exceeds \code{mer}. This holds
    under \code{pointit=TRUE} and \code{subset=n}. For \code{subset<n},
    the ratio is adjusted, but for small \code{subset}, the results
    may extremely vary and have to be taken with care.}
  \item{distcut}{optional value between 0 and 1. A similarity
    measure between FPCs, given in Hennig (2002), and the corresponding
    Single Linkage groups of FPCs with similarity larger
    than \code{distcut} are computed.
    A single representative FPC is selected for each group.}
  \item{maxit}{optional integer. Maximum number of iterations
    per algorithm run (usually an FPC is found much earlier).}
  \item{iter}{positive number. Algorithm stops when difference between
    subsequent weight vectors is smaller than \code{iter}. Only needed
    for \code{method="fuzzy"}.}
  \item{init.group}{optional list of logical vectors of length
    \code{n}.
    Every vector indicates a starting configuration for the fixed
    point algorithm. This can be used for datasets with high
    dimension, where the vectors of \code{init.group} indicate cluster
    candidates found by graphical inspection or background
    knowledge, as in Hennig and Christlieb (2002).}
  \item{ind.storage}{optional logical. If \code{TRUE},
    then all indicator
    vectors of found FPCs are given in the value of \code{fixmahal}.
    May need lots of memory, but is a bit faster.}
  \item{countmode}{optional positive integer. Every \code{countmode}
    algorithm runs \code{fixmahal} shows a message.}
  \item{plot}{optional string. If \code{"start"}, you get a scatterplot
    of the first two variables to highlight the initial configuration,
    \code{"iteration"} generates such a plot at each iteration,
    \code{"both"} does both (this may be very time consuming).
    The default is \code{"none"}.}
  \item{object}{object of class \code{mfpc}, output of \code{fixmahal}.}
  \item{x}{object of class \code{mfpc}, output of \code{fixmahal}.}
  \item{maxnc}{positive integer. Maximum number of FPCs
    to be reported.}
  \item{no}{positive integer. Number of the representative FPC to
    be plotted.}
  \item{bw}{optional logical. If \code{TRUE}, plot is black/white,
    FPC is
    indicated by different symbol. Else FPC is indicated red.}
  \item{main}{plot title.}
  \item{xlab}{label for x-axis. If \code{NULL}, a default text is used.}
  \item{ylab}{label for y-axis. If \code{NULL}, a default text is used.}
  \item{pch}{plotting symbol, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{col}{plotting color, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{gv}{logical vector (or, with \code{method="fuzzy"},
    vector of weights between 0 and 1) of length \code{n}.
    Indicates the initial
    configuration for the fixed point algorithm.} 
  \item{...}{additional parameters to be passed to \code{plot}
    (no effects elsewhere).}
}
\details{A (crisp) Mahalanobis FPC is a data subset 
  that reproduces itself under the following operation: \cr
  Compute mean and covariance matrix estimator for the data
  subset, and compute all points of the dataset for which the squared
  Mahalanobis distance is smaller than \code{ca}.\cr
  Fixed points of this operation can be considered as clusters,
  because they contain only
  non-outliers (as defined by the above mentioned procedure) and all other
  points are outliers w.r.t. the subset. \cr
  The current default is to compute fuzzy Mahalanobis FPCs, where the
  points in the subset have a membership weight between 0 and 1 and give
  rise to weighted means and covariance matrices.
  The new weights are then obtained by computing the weight function
  \code{\link{wfu}} of the squared Mahalanobis distances, i.e.,
  full weight for squared distances smaller than \code{ca}, zero
  weight for squared distances larger than \code{ca2} and 
  decreasing weights (linear function of squared distances)
  in between.\cr
  A fixed point algorithm is started from the whole dataset,
  algorithms are started from the subsets specified in
  \code{init.group}, and further algorithms are started from further
  initial configurations as explained under \code{subset} and in the
  function \code{\link{mahalconf}}. \cr
  Usually some of the FPCs are unstable, and more than one FPC may
  correspond to the same significant pattern in the data. Therefore the
  number of FPCs is reduced: A similarity matrix is computed
  between FPCs. Similarity between sets is defined as the ratio between
  2 times size of
  intersection and the sum of sizes of both sets. The Single Linkage
  clusters (groups)
  of level \code{distcut} are computed, i.e. the connectivity
  components of the graph where edges are drawn between FPCs with
  similarity larger than \code{distcut}. Groups of FPCs whose members
  are found often enough (cf. parameter \code{mer}) are considered as
  stable enough. A representative FPC is
  chosen for every Single Linkage cluster of FPCs according to the
  maximum expectation ratio \code{ser}. \code{ser} is the ratio between
  the number of findings of an FPC and the number of points
  of an FPC, adjusted suitably if \code{subset<n}.
  Usually only the representative FPCs of stable groups
  are of interest. \cr
  Default tuning constants are taken from Hennig (2005).\cr
  Generally, the default settings are recommended for
  \code{fixmahal}. For large datasets, the use of
  \code{init.group} together with \code{pointit=FALSE}
  is useful. Occasionally, \code{mnc} and \code{startn} may be chosen
  smaller than the default,
  if smaller clusters are of interest, but this may lead to too many
  clusters. Decrease of
  \code{ca} will often lead to too many clusters, even for homogeneous
  data. Increase of \code{ca} will produce only very strongly
  separated clusters. Both may be of interest occasionally.\cr
  Singular covariance matrices during the iterations are handled by
  \code{\link{solvecov}}.

  \code{summary.mfpc} gives a summary about the representative FPCs of
  stable groups.

  \code{plot.mfpc} is a plot method for the representative FPC of stable
  group no. \code{no}. It produces a scatterplot, where
  the points belonging to the FPC are highlighted, the mean is and
  for \code{p<3} also the region of the FPC is shown. For \code{p>=3},
  the optimal separating projection computed by \code{\link{batcoord}}
  is shown.  

  \code{fpclusters.mfpc} produces a list of indicator vectors for the
  representative FPCs of stable groups.

  \code{fpmi} is called by \code{fixmahal} for a single fixed point
  algorithm and will usually not be executed alone.
}
\value{
  \code{fixmahal} returns an object of class \code{mfpc}. This is a list
  containing the components \code{nc, g, means, covs, nfound, er, tsc,
    ncoll, skc, grto, imatrix, smatrix, stn, stfound, ser, sfpc, ssig,
    sto, struc, n, p, method, cgen, ca, ca2, cvec, calpha, pointit,
    subset, mnc, startn, mer, distcut}.

  \code{summary.mfpc} returns an object of class \code{summary.mfpc}.
  This is a list containing the components \code{means, covs, stn,
    stfound, sn, ser, tskip, skc, tsc, sim, ca, ca2, calpha, mer, method,
    cgen, pointit}.

  \code{fpclusters.mfpc} returns a list of indicator vectors for the
  representative FPCs of stable groups.

  
  \code{fpmi} returns a list with the components \code{mg, covg, md,
    gv, coll, method, ca}.
  
  \item{nc}{integer. Number of FPCs.}
  \item{g}{list of logical vectors. Indicator vectors of FPCs. \code{FALSE}
    if \code{ind.storage=FALSE}.}
  \item{means}{list of numerical vectors. Means of FPCs. In
    \code{summary.mfpc}, only for representative
    FPCs of stable groups and sorted according to
    \code{ser}.}
  \item{covs}{list of numerical matrices. Covariance matrices of FPCs. In
    \code{summary.mfpc}, only for representative
    FPCs of stable groups and sorted according to
    \code{ser}.}
  \item{nfound}{vector of integers. Number of findings for the FPCs.}
  \item{er}{numerical vector. Ratio of number of findings of FPCs to their
    size. Under \code{pointit=TRUE},
    this can be taken as a measure of stability of FPCs.}
  \item{tsc}{integer. Number of algorithm runs leading to too small or
    too seldom found FPCs.}
  \item{ncoll}{integer. Number of algorithm runs where collinear
    covariance matrices occurred.}
  \item{skc}{integer. Number of skipped clusters.}
  \item{grto}{vector of integers. Numbers of FPCs to which algorithm
    runs led, which were started by \code{init.group}.}
  \item{imatrix}{vector of integers. Size of intersection between
    FPCs. See \code{\link{sseg}}.}
  \item{smatrix}{numerical vector. Similarities between
    FPCs. See \code{\link{sseg}}.}
  \item{stn}{integer. Number of representative FPCs of stable groups.
    In \code{summary.mfpc}, sorted according to \code{ser}.}
  \item{stfound}{vector of integers. Number of findings of members of
    all groups of FPCs. In \code{summary.mfpc}, sorted according to
    \code{ser}.}
  \item{ser}{numerical vector. Ratio of number of findings of groups of
    FPCs to their size. Under \code{pointit=TRUE},
    this can be taken as a measure of stability of FPCs. In
    \code{summary.mfpc}, sorted from largest to smallest.}
  \item{sfpc}{vector of integers.
    Numbers of representative FPCs of all groups.}
  \item{ssig}{vector of integers of length \code{stn}.
    Numbers of representative FPCs of the stable groups.}
  \item{sto}{vector of integers. Numbers of groups ordered
    according to largest \code{ser}.}
  \item{struc}{vector of integers. Number of group an FPC belongs to.}
  \item{n}{see arguments.}
  \item{p}{see arguments.}
  \item{method}{see arguments.}
  \item{cgen}{see arguments.}
  \item{ca}{see arguments, if \code{cgen} has been \code{"fixed"}. Else
    numerical vector of length \code{nc} (see below), giving the
    final values of \code{ca} for all FPC. In \code{fpmi}, tuning
    constant for the iterated FPC.}
  \item{ca2}{see arguments.}
  \item{cvec}{numerical vector of length \code{n} for
    \code{cgen="auto"}. The values for the
    tuning constant \code{ca} corresponding to the cluster sizes from
    \code{1} to \code{n}.}
  \item{calpha}{see arguments.}
  \item{pointit}{see arguments.}
  \item{subset}{see arguments.}
  \item{mnc}{see arguments.}
  \item{startn}{see arguments.}
  \item{mer}{see arguments.}
  \item{distcut}{see arguments.}

  \item{sn}{vector of integers. Number of points of representative
    FPCs.}
  \item{tskip}{integer. Number of algorithm runs leading to skipped FPCs.}
  \item{sim}{vector of integers. Size of intersections between
    representative FPCs of stable groups. See \code{\link{sseg}}.}

  \item{mg}{mean vector.}
  \item{covg}{covariance matrix.}
  \item{md}{Mahalanobis distances.}
  \item{gv}{logical (numerical, respectively, if \code{method="fuzzy"})
    indicator vector of iterated FPC.}
  \item{coll}{logical. \code{TRUE} means that singular covariance
    matrices occurred during the iterations.}
}
\references{
  Hennig, C. (2002) Fixed point clusters for linear regression:
  computation and comparison, \emph{Journal of
    Classification} 19, 249-276.

  Hennig, C. (2005) Fuzzy and Crisp Mahalanobis Fixed Point Clusters,
  in Baier, D., Decker, R., and Schmidt-Thieme, L. (eds.):
  \emph{Data Analysis and Decision Support}. Springer, Heidelberg,
  47-56, \url{http://www.homepages.ucl.ac.uk/~ucakche/papers/fuzzyfix.ps}
  

  Hennig, C. and Christlieb, N. (2002) Validating visual clusters in
  large datasets: Fixed point clusters of spectral features,
  \emph{Computational Statistics and Data Analysis} 40, 723-739.

}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\seealso{
\code{\link{fixreg}} for linear regression fixed point clusters.

\code{\link{mahalconf}}, \code{\link{wfu}}, \code{\link{cmahal}}
for computation of initial configurations, weights, tuning constants.  

\code{\link{sseg}} for indexing the similarity/intersection vectors
computed by \code{fixmahal}.

\code{\link{batcoord}}, \code{\link{cov.rob}}, \code{\link{solvecov}},
\code{\link{cov.wml}}, \code{\link{plotcluster}}
for computation of projections, (inverted)
covariance matrices, plotting.

\code{\link{rFace}} for generation of example data, see below.
}
\examples{
  set.seed(20000)
  face <- rFace(400,dMoNo=2,dNoEy=0, p=3)
  # The first example uses grouping information via init.group.
  initg <- list()
  grface <- as.integer(attr(face,"grouping"))
  for (i in 1:5) initg[[i]] <- (grface==i)
  ff0 <- fixmahal(face, pointit=FALSE, init.group=initg)
  summary(ff0)
  cff0 <- fpclusters(ff0)
  plot(face, col=1+cff0[[1]])
  plot(face, col=1+cff0[[4]]) # Why does this come out as a cluster? 
  plot(ff0, face, 4) # A bit clearer...
  # Without grouping information, examples need more time:
  # ff1 <- fixmahal(face)
  # summary(ff1)
  # cff1 <- fpclusters(ff1)
  # plot(face, col=1+cff1[[1]])
  # plot(face, col=1+cff1[[6]]) # Why does this come out as a cluster? 
  # plot(ff1, face, 6) # A bit clearer...
  # ff2 <- fixmahal(face,method="ml")
  # summary(ff2)
  # ff3 <- fixmahal(face,method="ml",calpha=0.95,subset=50)
  # summary(ff3)
  ## ...fast, but lots of clusters. mer=0.3 may be useful here.
  # set.seed(3000)
  # face2 <- rFace(400,dMoNo=2,dNoEy=0)
  # ff5 <- fixmahal(face2)
  # summary(ff5)
  ## misses right eye of face data; with p=6,
  ## initial configurations are too large for 40 point clusters 
  # ff6 <- fixmahal(face2, startn=30)
  # summary(ff6)
  # cff6 <- fpclusters(ff6)
  # plot(face2, col=1+cff6[[3]])
  # plot(ff6, face2, 3)
  # x <- c(1,2,3,6,6,7,8,120)
  # ff8 <- fixmahal(x)
  # summary(ff8)
  # ...dataset a bit too small for the defaults...
  # ff9 <- fixmahal(x, mnc=3, startn=3)
  # summary(ff9)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}
\keyword{robust}% __ONLY ONE__ keyword per line


