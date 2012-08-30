\name{ridgeline.diagnosis}
\alias{ridgeline.diagnosis}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ridgeline plots, ratios and unimodality}
\description{
  Computes ridgeline ratios and unimodality checks for pairs of components
  given the parameters of a Gaussian mixture. Produces ridgeline plots.
}
\usage{
  ridgeline.diagnosis (propvector,muarray,Sigmaarray,
                                k=length(propvector),
                                ipairs="all", compute.ratio=TRUE,by=0.001,
                                ratiocutoff=NULL,ridgelineplot="matrix")

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{propvector}{vector of component proportions. Length must be
    number of components, and must sum up to 1.}
  \item{muarray}{matrix of component means (different components are in
    different columns).}
  \item{Sigmaarray}{three dimensional array with component covariance
    matrices (the third dimension refers to components).}
  \item{k}{integer. Number of components.}
  \item{ipairs}{\code{"all"} or list of vectors of two integers. If
    \code{ipairs="all"}, computations are carried out for all pairs of
    components. Otherwise, ipairs gives the pairs of components for
    which computations are carried out.}
  \item{compute.ratio}{logical. If \code{TRUE}, a matrix of ridgeline
    ratios is computed, see Hennig (2010a).}
  \item{by}{real between 0 and 1. Interval width for density computation
    along the ridgeline.}
  \item{ratiocutoff}{real between 0 and 1. If not \code{NULL}, the
    \code{connection.matrix} (see below) is computed by checking whether
    ridgeline ratios between components are below \code{ratiocutoff}.}
  \item{ridgelineplot}{one of \code{"none"}, \code{"matrix"},
    \code{"pairwise"}. If \code{"matrix"}, a matrix of pairwise
    ridgeline plots (see Hennig 2010b) will be plotted. If
    \code{"pairwise"}, pairwise ridgeline plots are plotted (you may
    want to set \code{par(ask=TRUE)} to see them all). No plotting if
    \code{"none"}.}
}

\value{
  A list with components
\item{merged.clusters}{vector of integers, stating for every mixture
  component the number of the cluster of components that would be merged
  by merging connectivity components of the graph specified by
  \code{connection.matrix}.}
\item{connection.matrix}{zero-one matrix, in which a one means that the
  mixture of the corresponding pair of components of the original
  mixture is either unimodel (if \code{ratiocutoff=NULL}) or that their
  ridgeline ratio is above \code{ratiocutoff}. If \code{ipairs!="all"},
  ignored pairs always have 0 in this matrix, same for
  \code{ratio.matrix}.}
\item{ratio.matrix}{matrix with entries between 0 und 1, giving the
  ridgeline ratio, which is the density minimum of the mixture of the
  corresponding pair of components along the ridgeline divided by the
  minimum of the two maxima closest to the beginning and the end of the
  ridgeline.}
}

\references{
  Hennig, C. (2010a) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
  
  Hennig, C. (2010b) Ridgeline plot and clusterwise stability as tools
  for merging Gaussian mixture components. To appear in
  \emph{Classification as a Tool for Research}, Proceedings of IFCS
  2009.
  
  Ray, S. and Lindsay, B. G. (2005) The Topography of Multivariate 
Normal Mixtures, \emph{Annals of Statistics}, 33, 2042-2065.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{ridgeline}}, \code{\link{dridgeline}},
  \code{\link{piridge}}, \code{\link{piridge.zeroes}}
}

\examples{
  muarray <- cbind(c(0,0),c(0,0.1),c(10,10))
  sigmaarray <- array(c(diag(2),diag(2),diag(2)),dim=c(2,2,3))
  rd <-
  ridgeline.diagnosis(c(0.5,0.3,0.2),muarray,sigmaarray,ridgelineplot="matrix",by=0.1)
  # Much slower but more precise with default by=0.001.
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



