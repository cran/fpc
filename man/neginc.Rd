\name{neginc}
\alias{neginc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Neg-entropy normality index for cluster validation}
\description{
  Cluster validity index based on the neg-entropy distances of
  within-cluster distributions to normal distribution, see
  Lago-Fernandez and Corbacho (2010). 
}
\usage{
neginc(x,clustering)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{something that can be coerced into a numerical
  matrix. Euclidean dataset.}
  \item{clustering}{vector of integers with length \code{=nrow(x)};
  indicating the cluster for each observation.}
}

\value{
  Index value, see
  Lago-Fernandez and Corbacho (2010). The lower (i.e., the more
  negative) the better. 
}

\references{
  Lago-Fernandez, L. F. and Corbacho, F. (2010) Normality-based
  validation for crisp clustering. \emph{Pattern Recognition} 43, 782-795.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
  options(digits=3)
  iriss <- as.matrix(iris[c(1:10,51:55,101:105),-5])
  irisc <- as.numeric(iris[c(1:10,51:55,101:105),5])
  neginc(iriss,irisc)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




