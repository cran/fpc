\name{plotcluster}
\alias{plotcluster}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Discriminant projection plot.}
\description{
  Plots to distinguish given classes by ten available projection
  methods. Includes classical discriminant
  coordinates, methods to project differences in
  mean and covariance structure, asymmetric methods (separation of a
  homogeneous class from a heterogeneous one), local neighborhood-based
  methods and methods based on robust covariance matrices.
  One-dimensional data is plotted against the cluster number.
}
\usage{
plotcluster(x, clvecd, clnum=NULL,
            method=ifelse(is.null(clnum),"dc","awc"),
            bw=FALSE,
            ignorepoints=FALSE, ignorenum=0, pointsbyclvecd=TRUE,
            xlab=NULL, ylab=NULL,
            pch=NULL, col=NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{vector of class numbers which can be coerced into
    integers; length must equal
    \code{nrow(xd)}.}
  \item{method}{one of
    \describe{
      \item{"dc"}{usual discriminant coordinates, see \code{\link{discrcoord}},}
      \item{"bc"}{Bhattacharyya coordinates, first coordinate showing
	mean differences, second showing covariance matrix differences,
	see \code{\link{batcoord}},}
      \item{"vbc"}{variance dominated Bhattacharyya coordinates,
	see \code{\link{batcoord}},}
      \item{"mvdc"}{added mean and variance differences optimizing
	coordinates, see \code{\link{mvdcoord}},}
      \item{"adc"}{asymmetric discriminant coordinates, see
	\code{\link{adcoord}},}
      \item{"awc"}{asymmetric discriminant coordinates with weighted
	observations, see \code{\link{awcoord}},}
      \item{"arc"}{asymmetric discriminant coordinates with weighted
	observations and robust MCD-covariance matrix,
	see \code{\link{awcoord}},}
      \item{"nc"}{neighborhood based coordinates,
	see \code{\link{ncoord}},}
      \item{"wnc"}{neighborhood based coordinates with weighted neighborhoods,
	see \code{\link{ncoord}},}
      \item{"anc"}{asymmetric neighborhood based coordinates,
	see \code{\link{ancoord}}.}
    }
    Note that "bc", "vbc", "adc", "awc", "arc" and "anc" assume that
    there are only two classes.}
  \item{clnum}{integer. Number of the class which is attempted to plot
    homogeneously by "asymmetric methods", which are the methods
    assuming that there are only two classes, as indicated above.
    \code{clnum} is ignored for methods "dc" and "nc".}
  \item{bw}{logical. If \code{TRUE}, the classes are distinguished by
    symbols, and the default color is black/white.
    If \code{FALSE}, the classes are distinguished by
    colors, and the default symbol is \code{pch=1}.}
  \item{ignorepoints}{logical. If \code{TRUE}, points with label
    \code{ignorenum} in \code{clvecd} are ignored in the computation for
    \code{method} and are only projected afterwards onto the resulting
    units. If \code{pch=NULL}, the plot symbol for these points is "N".}
  \item{ignorenum}{one of the potential values of the components of
    \code{clvecd}. Only has effect if \code{ignorepoints=TRUE}, see above.}
  \item{pointsbyclvecd}{logical. If \code{TRUE} and \code{pch=NULL}
    and/or \code{col=NULL}, some hopefully suitable
    plot symbols (numbers and letters) and colors are chosen to
    distinguish the values of \code{clvecd}, starting with "1"/"black"
    for the cluster with the smallest \code{clvecd}-code (note that
    colors for clusters with numbers larger than minimum number
    \code{+3} are drawn at random from all available colors).
    \code{FALSE} produces
    potentially less reasonable (but nonrandom) standard colors and symbols if
    \code{method} is "dc" or "nc", and will only distinguish whether
    \code{clvecd=clnum} or not for the other methods.} 
  \item{xlab}{label for x-axis. If \code{NULL}, a default text is used.}
  \item{ylab}{label for y-axis. If \code{NULL}, a default text is used.}
  \item{pch}{plotting symbol, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{col}{plotting color, see \code{\link{par}}.
    If \code{NULL}, the default is used.}
  \item{...}{additional parameters passed to \code{plot} or the
    projection methods.}
}
% \details{

% }
\note{
  For some of the asymmetric methods, the area in the plot
  occupied by the "homogeneous class" (see \code{clnum} above) may be
  very small, and it may make sense to run \code{plotcluster} a second
  time specifying plot parameters \code{xlim} and \code{ylim} in a
  suitable way. It often makes sense to magnify the plot region
  containing the homogeneous class in this way
  so that its separation from the rest can be
  seen more clearly.
  }
  

\references{
  Hennig, C. (2004) Asymmetric linear dimension reduction for classification.
  Journal of Computational and Graphical Statistics 13, 930-945 .
  
  Hennig, C. (2005)  A method for visual cluster validation.  In:
  Weihs, C. and Gaul, W. (eds.): Classification - The Ubiquitous
  Challenge. Springer, Heidelberg 2005, 153-160.
  
  Seber, G. A. F. (1984). \emph{Multivariate Observations}. New York: Wiley.

  Fukunaga (1990). \emph{Introduction to Statistical Pattern
  Recognition} (2nd ed.). Boston: Academic Press.
}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}}

\seealso{
  \code{\link{discrcoord}}, \code{\link{batcoord}},
  \code{\link{mvdcoord}}, \code{\link{adcoord}},
  \code{\link{awcoord}}, \code{\link{ncoord}},
  \code{\link{ancoord}}.

  \code{\link{discrproj}} is an interface to all these projection methods.

  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
set.seed(4634)
face <- rFace(300,dMoNo=2,dNoEy=0)
grface <- as.integer(attr(face,"grouping"))
plotcluster(face,grface)
plotcluster(face,grface==1)
plotcluster(face,grface, clnum=1, method="vbc")
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



