\name{cluster.varstats}
\alias{cluster.varstats}
\alias{print.varwisetables}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Variablewise statistics for clusters}
\description{
  This function gives some helpful variable-wise information for cluster
  interpretation, given a clustering and a data set. The output object
  contains some tables. For categorical variables, tables compare
  clusterwise distributions with overall distributions. Continuous
  variables are categorised for this.

  If desired, tables, histograms, some standard statistics of
  continuous variables and validation plots as available through
  \code{\link{discrproj}} (Hennig 2004) are given out on the fly. 
}
\usage{
cluster.varstats(clustering,vardata,contdata=vardata,
                             clusterwise=TRUE,
                            tablevar=NULL,catvar=NULL,
                             quantvar=NULL, catvarcats=10,
                             proportions=FALSE,
                            projmethod="none",minsize=ncol(contdata)+2,
                          ask=TRUE,rangefactor=1)

\method{print}{varwisetables}(x,digits=3,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{clustering}{vector of integers. Clustering (needs to be in
    standard coding, 1,2,...).}
  \item{vardata}{data matrix or data frame of which variables are
    summarised.}
  \item{contdata}{variable matrix or data frame, normally all or some
    variables from \code{vardata}, on which cluster visualisation by
    projection methods is performed unless \code{projmethod="none"}. It
    should make sense to interpret these variables in a quantitative
    (interval-scaled) way.}
  \item{clusterwise}{logical. If \code{FALSE}, only the output tables
    are computed but no more detail and graphs are given on the fly.}
  \item{tablevar}{vector of integers. Numbers of variables treated as
    categorical (i.e., no histograms and statistics, just tables) if
    \code{clusterwise=TRUE}. Note
    that an error will be produced by factor type variables unless they
    are declared as categorical here.}
  \item{catvar}{vector of integers. Numbers of variables to be
    categorised by proportional quantiles for table computation.
    Recommended for all continuous variables.}
  \item{quantvar}{vector of integers. Variables for which means,
    standard deviations and quantiles should be given out if
    \code{clusterwise=TRUE}.}
  \item{catvarcats}{integer. Number of categories used for
    categorisation of variables specified in \code{quantvar}.}
  \item{proportions}{logical. If \code{TRUE}, output tables contain
    proportions, otherwise numbers of observations.}
  \item{projmethod}{one of \code{"none"}, \code{"dc"}, \code{"bc"},
    \code{"vbc"}, \code{"mvdc"}, \code{"adc"}, \code{"awc"} (recommended
    if not \code{"none"}), \code{"arc"}, \code{"nc"}, \code{"wnc"},
    \code{"anc"}. Cluster validation projection method introduced in
    Hennig (2004), passed on as \code{method} argument in
    \code{\link{discrproj}}.}
  \item{minsize}{integer. Projection is not carried out for clusters
    with fewer points than this. (If this is chosen smaller, it may lead
    to errors with some projection methods.)}
  \item{ask}{logical. If \code{TRUE}, \code{par(ask=TRUE)} is set in the
    beginning to prompt the user before plots and \code{par(ask=FALSE)}
    in the end.}
  \item{rangefactor}{numeric. Factor by which to multiply the range for
    projection plot ranges.}
  \item{x}{an object of class \code{"varwisetables"}, output object of
    \code{cluster.varstats}.}
  \item{digits}{integer. Number of digits after the decimal point to
    print out.}
  \item{...}{not used.}
}

\value{
  An object of class \code{"varwisetables"}, which is a
  list with a table for each variable, giving (categorised) marginal
  distributions by cluster.
}

\references{
  Hennig, C. (2004) Asymmetric linear dimension reduction for classification.
  Journal of Computational and Graphical Statistics 13, 930-945 .
  }

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\examples{
  set.seed(112233)
  options(digits=3)
  require(MASS)
  require(flexmix)
  data(Cars93)
  Cars934 <- Cars93[,c(3,5,8,10)]
  cc <-
    discrete.recode(Cars934,xvarsorted=FALSE,continuous=c(2,3),discrete=c(1,4))
  fcc <- flexmix(cc$data~1,k=2,
  model=lcmixed(continuous=2,discrete=2,ppdim=c(6,3),diagonal=TRUE))
  cv <-
    cluster.varstats(fcc@cluster,Cars934, contdata=Cars934[,c(2,3)],
    tablevar=c(1,4),catvar=c(2,3),quantvar=c(2,3),projmethod="awc",
    ask=FALSE)
  print(cv)
}

\keyword{cluster}% __ONLY ONE__ keyword per line
