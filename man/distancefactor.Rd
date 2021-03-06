\name{distancefactor}
\alias{distancefactor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Factor for dissimilarity of mixed type data}
\description{
  Computes a factor that can be used to standardise ordinal categorical
  variables and binary dummy variables coding categories of nominal scaled
  variables for Euclidean
  dissimilarity computation in mixed type data. See Hennig and Liao (2013).
}
\usage{
distancefactor(cat,n=NULL, catsizes=NULL,type="categorical",
               normfactor=2,qfactor=ifelse(type=="categorical",1/2,
                             1/(1+1/(cat-1))))

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{cat}{integer. Number of categories of the variable to be standardised.
    Note that for \code{type="categorical"} the number of categories of
    the original variable is required, although the
    \code{distancefactor} is used to standardise dummy
    variables for the categories.}
  \item{n}{integer. Number of data points.}
  \item{catsizes}{vector of integers giving numbers of observations per
    category. One of \code{n} and \code{catsizes} must be supplied. If
    \code{catsizes=NULL}, \code{rep(round(n/cat),cat)} is used (this may
    be appropriate as well if numbers of observations of categories are
    unequal, if the researcher decides that the dissimilarity measure
    should not be influenced by empirical category sizes.} 
   \item{type}{\code{"categorical"} if the factor is used for dummy
     variables belonging to a nominal variable, \code{"ordinal"} if the
     factor is used for an ordinal variable ind standard Likert coding.}
   \item{normfactor}{numeric. Factor on which standardisation is based.
     As a default, this is \code{E(X_1-X_2)^2=2} for independent unit
     variance variables.}
   \item{qfactor}{numeric. Factor q in Hennig and Liao (2013) to
     adjust for clumping effects due to discreteness.}
}

\value{
  A factor by which to multiply the variable in order to make it
  comparable to a unit variance continuous variable when aggregated in
  Euclidean fashion for dissimilarity computation, so that expected
  effective difference between two realisations of the variable equals
  \code{qfactor*normfactor}. 
}


\references{
  Hennig, C. and Liao, T. (2013) How to find an appropriate clustering
  for mixed-type variables with application to socio-economic
  stratification, \emph{Journal of the Royal Statistical Society, Series
  C Applied Statistics}, 62, 309-369.


  }

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{\code{\link{lcmixed}}, \code{\link[cluster]{pam}}}

\examples{
  set.seed(776655)
  d1 <- sample(1:5,20,replace=TRUE)
  d2 <- sample(1:4,20,replace=TRUE)
  ldata <- cbind(d1,d2)
  lc <- cat2bin(ldata,categorical=1)$data
  lc[,1:5] <- lc[,1:5]*distancefactor(5,20,type="categorical")
  lc[,6] <- lc[,6]*distancefactor(4,20,type="ordinal")
}

\keyword{cluster}% __ONLY ONE__ keyword per line
