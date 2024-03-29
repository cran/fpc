\name{randomclustersim}
\alias{randomclustersim}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Simulation of validity indexes based on random clusterings}
\description{
  For a given dataset this simulates random clusterings using
  \code{\link{stupidkcentroids}}, \code{\link{stupidknn}},
  \code{\link{stupidkfn}}, and \code{\link{stupidkaven}}. It then
  computes and stores a set of cluster validity indexes for every
  clustering. 
}
\usage{
  randomclustersim(datadist,datanp=NULL,npstats=FALSE,useboot=FALSE,
                      bootmethod="nselectboot",
                      bootruns=25, 
                      G,nnruns=100,kmruns=100,fnruns=100,avenruns=100,
                      nnk=4,dnnk=2,
                      pamcrit=TRUE, 
                      multicore=FALSE,cores=detectCores()-1,monitor=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{datadist}{distances on which validation-measures are based, \code{dist}
    object or distance matrix.}
  \item{datanp}{optional observations times variables data matrix, see
    \code{npstats}.} 
  \item{npstats}{logical. If \code{TRUE}, \code{\link{distrsimilarity}}
    is called and the two statistics computed there are added to the
    output. These are based on \code{datanp} and require \code{datanp}
    to be specified.}
  \item{useboot}{logical. If \code{TRUE}, a stability index (either
    \code{nselectboot} or \code{prediction.strength}) will be involved.}
  \item{bootmethod}{either \code{"nselectboot"} or
    \code{"prediction.strength"}; stability index to be used if
    \code{useboot=TRUE}.}
  \item{bootruns}{integer. Number of resampling runs. If
    \code{useboot=TRUE}, passed on as \code{B} to \code{nselectboot} or
    \code{M} to \code{prediction.strength}.}
  \item{G}{vector of integers. Numbers of clusters to consider.}
  \item{nnruns}{integer. Number of runs of \code{\link{stupidknn}}.}
  \item{kmruns}{integer. Number of runs of \code{\link{stupidkcentroids}}.}
  \item{fnruns}{integer. Number of runs of \code{\link{stupidkfn}}.}
  \item{avenruns}{integer. Number of runs of \code{\link{stupidkaven}}.}
  \item{nnk}{\code{nnk}-argument to be passed on to
    \code{\link{cqcluster.stats}}.} 
  \item{dnnk}{\code{nnk}-argument to be passed on to
    \code{\link{distrsimilarity}}.}
  \item{pamcrit}{\code{pamcrit}-argument to be passed on to
    \code{\link{cqcluster.stats}}.} 
  \item{multicore}{logical. If \code{TRUE}, parallel computing is used
    through the function \code{\link{mclapply}} from package
    \code{parallel}; read warnings there if you intend to use this; it
    won't work on Windows.}
  \item{cores}{integer. Number of cores for parallelisation.}
  \item{monitor}{logical. If \code{TRUE}, it will print some runtime
    information.} 
}

% \details{
% }
\value{
  List with components
  \item{nn}{list, indexed by number of clusters. Every entry is
    a data frame with \code{nnruns} observations for every simulation
    run of \code{\link{stupidknn}}. The variables of the data frame are
    \code{avewithin, mnnd,
      cvnnd, maxdiameter, widestgap, sindex, minsep, asw, dindex,
      denscut, highdgap, pearsongamma, withinss, entropy}, if
    \code{pamcrit=TRUE} also \code{pamc}, if \code{npstats=TRUE} also
    \code{kdnorm, kdunif}. All these are cluster validation indexes;
    documented as values of \code{\link{clustatsum}}.}
  \item{fn}{list, indexed by number of clusters. Every entry is
    a data frame with \code{fnruns} observations for every simulation
    run of \code{\link{stupidkfn}}. The variables of the data frame are
    \code{avewithin, mnnd,
      cvnnd, maxdiameter, widestgap, sindex, minsep, asw, dindex,
      denscut, highdgap, pearsongamma, withinss, entropy}, if
    \code{pamcrit=TRUE} also \code{pamc}, if \code{npstats=TRUE} also
    \code{kdnorm, kdunif}. All these are cluster validation indexes;
    documented as values of \code{\link{clustatsum}}.}
  \item{aven}{list, indexed by number of clusters. Every entry is
    a data frame with \code{avenruns} observations for every simulation
    run of \code{\link{stupidkaven}}. The variables of the data frame are
    \code{avewithin, mnnd,
      cvnnd, maxdiameter, widestgap, sindex, minsep, asw, dindex,
      denscut, highdgap, pearsongamma, withinss, entropy}, if
    \code{pamcrit=TRUE} also \code{pamc}, if \code{npstats=TRUE} also
    \code{kdnorm, kdunif}. All these are cluster validation indexes;
    documented as values of \code{\link{clustatsum}}.}
  \item{km}{list, indexed by number of clusters. Every entry is
    a data frame with \code{kmruns} observations for every simulation
    run of \code{\link{stupidkcentroids}}. The variables of the data
    frame are \code{avewithin, mnnd,
      cvnnd, maxdiameter, widestgap, sindex, minsep, asw, dindex,
      denscut, highdgap, pearsongamma, withinss, entropy}, if
    \code{pamcrit=TRUE} also \code{pamc}, if \code{npstats=TRUE} also
    \code{kdnorm, kdunif}. All these are cluster validation indexes;
    documented as values of \code{\link{clustatsum}}.}
  \item{nnruns}{number of involved runs of \code{\link{stupidknn}},}
  \item{fnruns}{number of involved runs of \code{\link{stupidkfn}},}
  \item{avenruns}{number of involved runs of \code{\link{stupidkaven}},}
  \item{kmruns}{number of involved runs of \code{\link{stupidkcentroids}},}
  \item{boot}{if \code{useboot=TRUE}, stability value; \code{stabk} for
    method \code{\link{nselectboot}}; \code{mean.pred} for method
    \code{\link{prediction.strength}}.}
  
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
  \code{\link{stupidkcentroids}}, \code{\link{stupidknn}}, \code{\link{stupidkfn}}, \code{\link{stupidkaven}}, \code{\link{clustatsum}}
}

\examples{
  set.seed(20000)
  options(digits=3)
  face <- rFace(10,dMoNo=2,dNoEy=0,p=2)
  rmx <- randomclustersim(dist(face),datanp=face,npstats=TRUE,G=2:3,
    nnruns=2,kmruns=2, fnruns=1,avenruns=1,nnk=2)
\dontrun{
  rmx$km # Produces slightly different but basically identical results on ATLAS
}
  rmx$aven
  rmx$fn
  rmx$nn
  
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line



