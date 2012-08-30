\name{con.comp}
\alias{con.comp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Connectivity components of an undirected graph}
\description{
  Computes the connectivity components of an undirected graph from a
  matrix giving the edges.
}
\usage{
con.comp(comat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{comat}{a symmetric logical or 0-1 matrix, where \code{comat[i,j]=TRUE}
    means that there is an edge between vertices \code{i} and
    \code{j}. The diagonal is ignored.}
}
\details{
  The "depth-first search" algorithm of Cormen, Leiserson and Rivest
  (1990, p. 477) is used.
}
\value{
  An integer vector, giving the number of the connectivity component for
  each vertice.
}
\references{
  Cormen, T. H., Leiserson, C. E. and Rivest, R. L. (1990), \emph{Introduction
    to Algorithms}, Cambridge: MIT Press.
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{hclust}}, \code{\link{cutree}} for cutted single linkage
  trees (often equivalent).
}

\examples{
  set.seed(1000)
  x <- rnorm(20)
  m <- matrix(0,nrow=20,ncol=20)
  for(i in 1:20)
    for(j in 1:20)
      m[i,j] <- abs(x[i]-x[j])
  d <- m<0.2
  cc <- con.comp(d)
  max(cc) # number of connectivity components
  plot(x,cc)
  # The same should be produced by
  # cutree(hclust(as.dist(m),method="single"),h=0.2).
}
\keyword{array}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line


