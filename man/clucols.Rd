\name{clucols}
\alias{clucols}
\alias{clugrey}
\alias{clusym}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sets of colours and symbols for cluster plotting}
\description{
  \code{clucols} gives out a vector of different random colours.
  \code{clugrey} gives out a vector of equidistant grey scales.
  \code{clusym} is a vector of different symbols starting from "1",
  "2",...
}
\usage{
  clucols(i, seed=NULL)
  clugrey(i,max=0.9)
  clusym
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{i}{integer. Length of output vector (number of clusters).}
  \item{seed}{integer. Random seed.}
  \item{max}{between 0 and 1. Maximum grey scale value, see
    \code{\link{grey}} (close to 1 is bright).}    
}

\value{
  \code{clucols} gives out a vector of different random colours.
  \code{clugrey} gives out a vector of equidistant grey scales.
  \code{clusym} is a vector of different characters starting from "1",
  "2",...
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\examples{
  set.seed(112233)
  data(Cars93)
  Cars934 <- Cars93[,c(3,5,8,10)]
  cc <-
    discrete.recode(Cars934,xvarsorted=FALSE,continuous=c(2,3),discrete=c(1,4))
  fcc <- flexmix(cc$data~1,k=3,
  model=lcmixed(continuous=2,discrete=2,ppdim=c(6,3),diagonal=TRUE))
  plot(Cars934[,c(2,3)],col=clucols(3)[fcc@cluster],pch=clusym[fcc@cluster])
}

\keyword{cluster}% __ONLY ONE__ keyword per line
