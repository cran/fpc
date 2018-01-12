\name{xtable}
\alias{xtable}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Partition crosstable with empty clusters}
\description{
  This produces a crosstable between two integer vectors (partitions) of
  the same length with a given maximum vector entry \code{k} so that the
  size of the table is \code{k*k} with zeroes for missing entries
  between 1 and \code{k} (the command \code{\link{table}} does pretty
  much the same thing but will leave out missing entries). 
}
\usage{
xtable(c1,c2,k)
}
\arguments{
  \item{c1}{vector of integers.}
  \item{c2}{vector of integers of same length as \code{c1}.}
  \item{k}{integer. Must be larger or equal to maximum entry in
    \code{c1} and \code{c2}.}
}

\value{
  A matrix of dimensions \code{c(k,k)}. Entry \code{[i,j]} gives the
  number of places in which \code{c1==i & c2==j}. 
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}

\seealso{
  \code{\link{table}}
}
\examples{
  c1 <- 1:3
  c2 <- c(1,1,2)
  xtable(c1,c2,3)
}
\keyword{array}% at least one, from doc/KEYWORDS



