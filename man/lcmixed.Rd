\name{lcmixed}
\alias{lcmixed}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{flexmix method for mixed Gaussian/multinomial mixtures}
\description{
  \code{lcmixed} is a method for the
  \code{\link[flexmix]{flexmix}}-function in package
  \code{flexmix}. It provides the necessary information to run an
  EM-algorithm for maximum likelihood estimation for a latent class
  mixture (clustering) model where some variables are continuous
  and modelled within the mixture components by Gaussian distributions
  and some variables are categorical and modelled within components by
  independent multinomial distributions. \code{lcmixed} can be called
  within \code{flexmix}. The function \code{\link{flexmixedruns}} is a wrapper
  function that can be run to apply \code{lcmixed}.

  Note that at least one categorical variable is needed, but it is
  possible to use data without continuous variable.

  There are further format restrictions to the data (see below in the
  documentation of \code{continuous} and \code{discrete}), which
  can be ignored when running \code{lcmixed} through
  \code{\link{flexmixedruns}}.
}
\usage{
lcmixed( formula = .~. , continuous, discrete, ppdim,
                     diagonal = TRUE, pred.ordinal=FALSE, printlik=FALSE )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a formula to specify response and explanatory
    variables. For \code{lcmixed} this always has the form \code{x~1},
    where \code{x} is a matricx or data frome of all variables to be
    involved, because regression and explanatory variables are not
    implemented.}
  \item{continuous}{number of continuous variables. Note that the
    continuous variables always need to be the first variables in the
    matrix or data frame.}
  \item{discrete}{number of categorical variables. Always the last
    variables in the matrix or data frame. Note that categorical
    variables always must be coded as integers 1,2,3, etc. without
    interruption.}
  \item{ppdim}{vector of integers specifying the number of (in the data)
    existing categories for each categorical variable.}
  \item{diagonal}{logical. If \code{TRUE}, Gaussian models are fitted
    restricted to diagonal covariance matrices. Otherwise, covariance
    matrices are unrestricted. \code{TRUE} is consistent with the
    "within class independence" assumption for the multinomial variables.}
  \item{pred.ordinal}{logical. If \code{FALSE}, the within-component
    predicted value for categorical variables is the probability mode,
    otherwise it is the mean of the standard (1,2,3,...) scores, which
    may be better for ordinal variables.}
  \item{printlik}{logical. If \code{TRUE}, the loglikelihood is printed
    out whenever computed.}
}

\details{
  The data need to be organised case-wise, i.e., if there are
  categorical variables only, and 15 cases with values c(1,1,2) on the
  3 variables, the data matrix needs 15 rows with values 1 1 2.
  
  General documentation on flexmix methods can be found in Chapter 4 of
  Friedrich Leisch's "FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R",
  \url{http://cran.r-project.org/web/packages/flexmix/vignettes/flexmix-intro.pdf}
}
 
\value{
  An object of class \code{FLXMC} (not documented; only used
  internally by \code{flexmix}).
}

\references{
  Hennig, C. and Liao, T. (2010) Comparing latent class and
  dissimilarity based clustering for mixed type variables with
  application to social stratification. Research report no. 308,
  Department of Statistical Science, UCL.
  \url{http://www.ucl.ac.uk/Stats/research/reports/psfiles/rr308.pdf}
  }

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{
  \code{\link{flexmixedruns}}, \code{\link[flexmix]{flexmix}},
  \code{\link[flexmix]{flexmix-class}},
  \code{\link{discrete.recode}}, which recodes a dataset into the format
  required by \code{lcmixed}
}

\examples{
  set.seed(112233)
  data(Cars93)
  Cars934 <- Cars93[,c(3,5,8,10)]
  cc <-
  discrete.recode(Cars934,xvarsorted=FALSE,continuous=c(2,3),discrete=c(1,4))
  fcc <- flexmix(cc$data~1,k=2,
  model=lcmixed(continuous=2,discrete=2,ppdim=c(6,3),diagonal=TRUE))
  summary(fcc)
}
\keyword{cluster}% __ONLY ONE__ keyword per line
