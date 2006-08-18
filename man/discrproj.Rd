\name{discrproj}
\alias{discrproj}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Linear dimension reduction for classification}
\description{
   An interface for ten methods of linear dimension reduction in order
   to separate the groups optimally in the projected data. Includes
   classical discriminant coordinates, methods to project differences in
   mean and covariance structure, asymmetric methods (separation of a
   homogeneous class from a heterogeneous one), local neighborhood-based
   methods and methods based on robust covariance matrices.
}
\usage{
 discrproj(x, clvecd, method="awc", clnum=1, ...)
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
      \item{"mvdc"}{added meana and variance differences optimizing
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
    assuming that there are only two classes, as indicated above.} 
  \item{...}{additional parameters passed to the
    projection methods.}
}
% \details{
% }

\value{
  \code{discrproj} returns the output of the chosen projection method,
  which is a list with at least the components \code{ev, units, proj}.
  For detailed informations see the help pages of the projection methods.
  \item{ev}{eigenvalues in descending order, usually indicating portion
    of information in the corresponding direction.}
  \item{units}{columns are coordinates of projection basis vectors.
  New points \code{x} can be projected onto the projection basis vectors
  by \code{x \%*\% units}}
  \item{proj}{projections of \code{xd} onto \code{units}.} 
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
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}}

\seealso{
  \code{\link{discrcoord}}, \code{\link{batcoord}},
  \code{\link{mvdcoord}}, \code{\link{adcoord}},
  \code{\link{awcoord}}, \code{\link{ncoord}},
  \code{\link{ancoord}}.

  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
set.seed(4634)
face <- rFace(300,dMoNo=2,dNoEy=0,p=3)
grface <- as.integer(attr(face,"grouping"))
discrproj(face,grface, method="nc")$units
discrproj(face,grface, method="wnc")$units
discrproj(face,grface==1, method="arc")$units
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line



