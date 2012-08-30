\name{mergeparameters}
\alias{mergeparameters}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{New parameters from merging two Gaussian mixture components}
\description{
  Re-computes pointwise posterior probabilities, mean and covariance
  matrix for a mixture component obtained by merging two mixture
  components in a Gaussian mixture.
}
\usage{
  mergeparameters(xdata, j1, j2, probs, muarray,Sigmaarray, z)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xdata}{data (something that can be coerced into a matrix).}
  \item{j1}{integer. Number of first mixture component to be merged.}
  \item{j2}{integer. Number of second mixture component to be merged.}
  \item{probs}{vector of component proportions (for all components;
    should sum up to one).}
  \item{muarray}{matrix of component means (rows).}
  \item{Sigmaarray}{array of component covariance matrices (third
    dimension refers to component number).}
  \item{z}{matrix of observation- (row-)wise posterior probabilities of
    belonging to the components (columns).} 
}

\value{
  List with components
  \item{probs}{see above; sum of probabilities for original components
  \code{j1} and \code{j2} is now \code{probs[j1]}. Note that generally,
  also for the further components, values for the merged component are
  in place \code{j1} and values in place \code{j2} are not changed. This
  means that in order to have only the information for the new mixture
  after merging, the entries in places \code{j2} need to be suppressed.}  
  \item{muarray}{see above; weighted mean of means of component
  \code{j1} and \code{j2} is now in place \code{j1}.}
  \item{Sigmaarray}{see above; weighted covariance matrix handled as
  above.} 
  \item{z}{see above; original entries for columns \code{j1} and
  \code{j2} are summed up and now in column \code{j1}.}
}

\references{
  Hennig, C. (2010) Methods for merging Gaussian mixture components,
  \emph{Advances in Data Analysis and Classification}, 4, 3-34.
}
\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}
\examples{
  set.seed(98765)
  iriss <- iris[sample(150,20),-5]
  irisBIC <- mclustBIC(iriss)
  siris <- summary(irisBIC,iriss)
  probs <- siris$parameters$pro
  muarray <- siris$parameters$mean
  Sigmaarray <- siris$parameters$variance$sigma
  z <- siris$z
  mergeparameters(iriss,1,2,probs,muarray,Sigmaarray,z)
}
\keyword{multivariate}
\keyword{cluster}


