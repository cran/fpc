\name{batcoord}
\alias{batcoord}
\alias{batvarcoord}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Bhattacharyya discriminant projection}
\description{
  Computes Bhattacharyya discriminant projection coordinates
  as described in Fukunaga (1990), p. 455 ff.
}
\usage{
batcoord(xd, clvecd, clnum=1, dom="mean")
batvarcoord(xd, clvecd, clnum=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xd}{the data matrix; a numerical object which can be coerced
    to a matrix.}
  \item{clvecd}{integer or logical vector of class numbers; length must equal
    \code{nrow(xd)}.}
  \item{clnum}{integer, one of the values of \code{clvecd}, if this is
    an integer vector. Bhattacharyya projections can only be computed if
    there are only two classes in the dataset. \code{clnum} is the number
    of one of the two classes. All the points indicated by other values
    of \code{clvecd} are interpreted as the second class.}
  \item{dom}{string. \code{dom="mean"} means that the discriminant
    coordinate for the group means is computed as the first projection
    direction by
    \code{\link{discrcoord}} (option \code{pool="equal"}; both classes
    have the same weight for computing the within-class covariance
    matrix). Then the data is projected into a subspace orthogonal
    (w.r.t. the within-class covariance) to the
    discriminant coordinate, and the projection coordinates to maximize
    the differences in variance are computed. \cr
    \code{dom="variance"} means that the projection coordinates
    maximizing the difference in variances are computed. Then they are
    ordered with respect to the Bhattacharyya distance, which takes also
    the mean differences into account. Both procedures are implemented
    as described in Fukunaga (1990).}
}
\details{
\code{batvarcoord} computes the optimal projection coordinates with
respect to the difference in variances. \code{batcoord} combines the
differences in mean and variance as explained for the argument \code{dom}.
}
\value{
  \code{batcoord} returns a list with the components \code{ev, rev,
    units, proj}. \code{batvarcoord} returns a list with the components
  \code{ev, rev, units, proj, W, S1, S2}. 
  \item{ev}{vector of eigenvalues. If \code{dom="mean"}, then first eigenvalue
    from \code{\link{discrcoord}}. Further eigenvalues are of
  \eqn{S_1^{-1}S_2}, where \eqn{S_i} is the covariance matrix of class
  i. For \code{batvarcoord} or
  if \code{dom="variance"}, all eigenvalues come from
  \eqn{S_1^{-1}S_2} and are ordered by \code{rev}.}
\item{rev}{for \code{batcoord}:
  vector of projected Bhattacharyya distances (Fukunaga
  (1990), p. 99). Determine quality of the projection coordinates.
  For \code{batvarcoord}: vector of amount of projected difference in
  variances.}
  \item{units}{columns are coordinates of projection basis vectors.
  New points \code{x} can be projected onto the projection basis vectors
  by \code{x \%*\% units}.}
  \item{proj}{projections of \code{xd} onto \code{units}.} 

  \item{W}{matrix \eqn{S_1^{-1}S_2}.}
  \item{S1}{covariance matrix of the first class.}    
  \item{S2}{covariance matrix of the second class.}    
}
\references{
  Fukunaga, K. (1990). \emph{Introduction to Statistical Pattern
    Recognition} (2nd ed.). Boston: Academic Press.
}
\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{
  \code{\link{plotcluster}} for straight forward discriminant plots.
  
  \code{\link{discrcoord}} for discriminant coordinates.

  \code{\link{rFace}} for generation of the example data used below.
}

\examples{
set.seed(4634)
face <- rFace(600,dMoNo=2,dNoEy=0)
grface <- as.integer(attr(face,"grouping"))
bcf2 <- batcoord(face,grface==2)
plot(bcf2$proj,col=1+(grface==2))
bcfv2 <- batcoord(face,grface==2,dom="variance")
plot(bcfv2$proj,col=1+(grface==2))
bcfvv2 <- batvarcoord(face,grface==2)
plot(bcfvv2$proj,col=1+(grface==2))
}
\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{classif}% __ONLY ONE__ keyword per line

