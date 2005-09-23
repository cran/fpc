\name{simmatrix}
\alias{simmatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Extracting intersections between clusters from fpc-object}
\description{
  Extracts the information about the size of the intersections
  between representative
  Fixed Point Clusters (FPCs) of stable groups from the output of
  the FPC-functions \code{\link{fixreg}} and \code{\link{fixmahal}}.
}
\usage{
simmatrix(fpcobj)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fpcobj}{an object of class \code{rfpc} or \code{mfpc}.}
}

\value{
  A non-negative real-valued vector giving the number of points in
  the intersections of the representative FPCs of stable groups.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\note{The intersection between representative FPCs no. \code{i} and
  \code{j} is at position \code{\link{sseg}(i,j)}.}

\seealso{
  \code{\link{fixmahal}},
  \code{\link{fixreg}},
  \code{\link{sseg}}
}

\examples{
set.seed(190000)
data(tonedata)
# Note: If you do not use the installed package, replace this by
# tonedata <- read.table("(path/)tonedata.txt", header=TRUE)
attach(tonedata)
tonefix <- fixreg(stretchratio,tuned,mtf=1,ir=20)
simmatrix(tonefix)[sseg(2,3)]
}
\keyword{utilities}% at least one, from doc/KEYWORDS

