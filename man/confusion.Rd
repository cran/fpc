\name{confusion}
\alias{confusion}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Misclassification probabilities in mixtures}
\description{
   Estimates a misclassification probability in a mixture distribution
   between two mixture components from estimated posterior probabilities
   regardless of component parameters, see Hennig (2010).
}
\usage{
confusion(z,pro,i,j,adjustprobs=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{z}{matrix of posterior probabilities for observations (rows) to
    belong to mixture components (columns), so entries need to sum up to
    1 for each row.}
  \item{pro}{vector of component proportions, need to sum up to 1.}
  \item{i}{integer. Component number.}
  \item{j}{integer. Component number.}
  \item{adjustprobs}{logical. If \code{TRUE}, probabilities are
    initially standardised so that those for components \code{i} and
    \code{j} add up to one (i.e., if they were the only components).}
}

\value{
  Estimated probability that an observation generated by component
  \code{j} is classified to component \code{i} by maximum a posteriori rule. 
}

\references{
  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
  }
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en/}
}
\examples{
  set.seed(12345)
  m <- rpois(20,lambda=5)
  dim(m) <- c(5,4)
  pro <- apply(m,2,sum)
  pro <- pro/sum(pro)
  m <- m/apply(m,1,sum)
  round(confusion(m,pro,1,2),digits=2)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{multivariate}



