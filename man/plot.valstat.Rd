\name{plot.valstat}
\alias{plot.valstat}
\alias{print.valstat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Simulation-standardised plot and print of cluster validation statistics}
\description{
  Visualisation and print function for cluster validation output
  compared to results
  on simulated random clusterings. The print method can also be used to
  compute and print an aggregated cluster validation index.

  Unlike for many other plot methods, the additional arguments
  of \code{plot.valstat} are essential. \code{print.valstat} should make
  good sense with the defaults, but for computing the aggregate index
   need to be set.
}
\usage{
\method{plot}{valstat}(x,simobject=NULL,statistic="sindex",
                            xlim=NULL,ylim=c(0,1),
                            nmethods=length(x)-5,
                            col=1:nmethods,cex=1,pch=c("c","f","a","n"),
                            simcol=rep(grey(0.7),4),
                         shift=c(-0.1,-1/3,1/3,0.1),include.othernc=NULL,...)


\method{print}{valstat}(x,statistics=x$statistics,
                          nmethods=length(x)-5,aggregate=FALSE,
                          weights=NULL,digits=2,
                          include.othernc=NULL,...)

			      
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{"valstat"}, such as sublists
    \code{stat, qstat, sstat} of \code{\link{clusterbenchstats}}-output.}
  \item{simobject}{list of simulation results as produced by
    \code{\link{randomclustersim}} and documented there; typically sublist
    \code{sim} of \code{\link{clusterbenchstats}}-output.}  
  \item{statistic}{one of \code{"avewithin","mnnd","variation",
      "diameter","gap","sindex","minsep","asw","dindex","denscut",
      "highdgap","pg","withinss","entropy","pamc","kdnorm","kdunif","dmode"};
    validation statistic to be plotted.}
  \item{xlim}{passed on to \code{plot}. Default is the range of all
      involved numbers of clusters, minimum minus 0.5 to maximum plus
      0.5.}
  \item{ylim}{passed on to \code{plot}.}    
  \item{nmethods}{integer. Number of clustering methods to involve
      (these are those from number 1 to \code{nmethods} specified in
      \code{x$name}).}
  \item{col}{colours used for the different clustering methods.}
  \item{cex}{passed on to \code{plot}.}
  \item{pch}{vector of symbols for random clustering results from
    \code{\link{stupidkcentroids}}, \code{\link{stupidkfn}},
    \code{\link{stupidkaven}}, \code{\link{stupidknn}}.
    To be passed on to \code{plot}.}
  \item{simcol}{vector of colours used for random clustering results in
    order \code{\link{stupidkcentroids}}, \code{\link{stupidkfn}},
    \code{\link{stupidkaven}}, \code{\link{stupidknn}}.}
  \item{shift}{numeric vector. Indicates the amount to which the results
    from  \code{\link{stupidkcentroids}}, \code{\link{stupidkfn}},
    \code{\link{stupidkaven}}, \code{\link{stupidknn}} are
    plotted to the right of their respective
    number of clusters (negative numbers plot to the left).}
  \item{include.othernc}{this indicates whether methods should be
    included that estimated their number of clusters themselves and gave
    a result outside the standard range as given by \code{x$minG}
    and \code{x$maxG}. If not \code{NULL}, this is a
    list of integer vectors of length 2. The first
    number is
    the number of the clustering method (the order is determined by
    argument \code{x$name}), the second number is the
    number of clusters for those methods that estimate the number of
    clusters themselves and estimated a number outside the standard
    range. Normally what will be used here, if not \code{NULL}, is the
    output parameter
    \code{cm$othernc} of \code{\link{clusterbenchstats}}, see also
    \code{\link{cluster.magazine}}.}
  \item{statistics}{vector of character strings specifying the
    validation statistics that will be included in the output (unless
    you want to restrict the output for some reason, the default should
    be fine.}
  \item{aggregate}{logical. If \code{TRUE}, an aggegate validation
    statistic will be computed as the weighted mean of the involved
    statistic. This requires \code{weights} to be set. In order for this
    to make sense, values of the validation statistics should be
    comparable, which is achieved by standardisation in
    \code{\link{clusterbenchstats}}. Accordingly, \code{x} should
    be the \code{qstat} or \code{sstat}-component of the
    \code{\link{clusterbenchstats}}-output rather than the
    \code{stat}-component.}
  \item{weights}{vector of numericals. Weights for computation of the
    aggregate statistic in case that \code{aggregate=TRUE}. The order of
    clustering methods corresponding to the weight vector is given by
    \code{x$name}.}
  \item{digits}{minimal number of significant digits, passed on to
    \code{\link{print.table}}.}
  \item{...}{no effect.}
}

\details{
  Whereas \code{print.valstat}, at least with \code{aggregate=TRUE}
  makes more sense for the \code{qstat} or \code{sstat}-component of the
  \code{\link{clusterbenchstats}}-output rather than the
  \code{stat}-component, \code{plot.valstat} should be run with the
  \code{stat}-component if \code{simobject} is specified, because the
  simulated cluster validity statistics are unstandardised and need to
  be compared with unstandardised values on the dataset of interest.

  \code{print.valstat} will print all values for all validation indexes
  and the aggregated index (in case of \code{aggregate=TRUE} and set
  \code{weights} will be printed last.
}
  

\value{
  \code{print.valstats} returns the results table as invisible object.
}

\references{
  Hennig, C. (2017) Cluster validation by measurement of clustering
  characteristics relevant to the user. In C. H. Skiadas (ed.)
  \emph{Proceedings of ASMDA 2017}, 501-520,
  \url{https://arxiv.org/abs/1703.09282}

  Akhanli, S. and Hennig, C. (2020) Calibrating and aggregating cluster
validity indexes for context-adapted comparison of clusterings. Accepted
for
publication by Statistics and Computing, \url{https://arxiv.org/abs/2002.01822}

}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\seealso{
  \code{\link{clusterbenchstats}}, \code{\link{valstat.object}}, 
  \code{\link{cluster.magazine}}
}
\examples{  
  set.seed(20000)
  options(digits=3)
  face <- rFace(10,dMoNo=2,dNoEy=0,p=2)
  clustermethod=c("kmeansCBI","hclustCBI","hclustCBI")
  clustermethodpars <- list()
  clustermethodpars[[2]] <- clustermethodpars[[3]] <- list()
  clustermethodpars[[2]]$method <- "ward.D2"
  clustermethodpars[[3]]$method <- "single"
  methodname <- c("kmeans","ward","single")
  cbs <-  clusterbenchstats(face,G=2:3,clustermethod=clustermethod,
    methodname=methodname,distmethod=rep(FALSE,3),
    clustermethodpars=clustermethodpars,nnruns=2,kmruns=2,fnruns=2,avenruns=2)
  plot(cbs$stat,cbs$sim)
  plot(cbs$stat,cbs$sim,statistic="dindex")
  plot(cbs$stat,cbs$sim,statistic="avewithin")
  print(cbs$sstat,aggregate=TRUE,weights=c(1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0))
# Some of the values are "NaN" because due to the low number of runs of
# the stupid clustering methods there is no variation. If this happens
# in a real application, nnruns etc. should be chosen higher than 2.
# Also useallg=TRUE in clusterbenchstats may help.  
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



