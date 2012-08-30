\name{flexmixedruns}
\alias{flexmixedruns}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Fitting mixed Gaussian/multinomial mixtures with flexmix}
\description{
  \code{flexmixedruns} fits a latent class
  mixture (clustering) model where some variables are continuous
  and modelled within the mixture components by Gaussian distributions
  and some variables are categorical and modelled within components by
  independent multinomial distributions. The fit is by maximum
  likelihood estimation computed with the EM-algorithm. The number of
  components can be estimated by the BIC.

  Note that at least one categorical variable is needed, but it is
  possible to use data without continuous variable.
}
\usage{
flexmixedruns(x,diagonal=TRUE,xvarsorted=TRUE,
                          continuous,discrete,ppdim=NULL,initial.cluster=NULL,
                          simruns=20,n.cluster=1:20,verbose=TRUE,recode=TRUE,
                          allout=TRUE,control=list(minprior=0.001),silent=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame. The data need to be organised
    case-wise, i.e., if there are categorical variables only, and 15
    cases with values c(1,1,2) on the 3 variables, the data matrix needs
    15 rows with values 1 1 2. (Categorical variables could take numbers
    or strings or anything that can be coerced to factor levels as values.)}
  \item{diagonal}{logical. If \code{TRUE}, Gaussian models are fitted
    restricted to diagonal covariance matrices. Otherwise, covariance
    matrices are unrestricted. \code{TRUE} is consistent with the
    "within class independence" assumption for the multinomial
    variables.}
  \item{xvarsorted}{logical. If \code{TRUE}, the continuous variables
    are assumed to be the first ones, and the categorical variables to
    be behind them.}
  \item{continuous}{vector of integers giving positions of the
    continuous variables. If \code{xvarsorted=TRUE}, a single integer,
    number of continuous variables.}
  \item{discrete}{vector of integers giving positions of the
    categorical variables. If \code{xvarsorted=TRUE}, a single integer,
    number of categorical variables.}
  \item{ppdim}{vector of integers specifying the number of (in the data)
    existing categories for each categorical variable. If
    \code{recode=TRUE}, this can be omitted and is computed
    automatically.}
  \item{initial.cluster}{this corresponds to the \code{cluster}
    parameter in \code{flexmix} and should only be specified if
    \code{simruns=1} and \code{n.cluster} is a single number.
    Either a matrix with \code{n.cluster} columns of initial cluster
    membership probabilities for each observation; or a factor or
    integer vector with the initial cluster assignments of
    observations at the start of the EM algorithm. Default is
    random assignment into \code{n.cluster} clusters.}
  \item{simruns}{integer. Number of starts of the EM algorithm with
    random initialisation in order to find a good global optimum.}
  \item{n.cluster}{vector of integers, numbers of components (the optimum
    one is found by minimising the BIC).}
  \item{verbose}{logical. If \code{TRUE}, some information about the
    different runs of the EM algorithm is given out.}
  \item{recode}{logical. If \code{TRUE}, the function
    \code{discrete.recode} is applied in order to recode categorical
    data so that the \code{lcmixed}-method can use it. Only set this
    to \code{FALSE} if your data already has that format (even it that
    case, \code{TRUE} doesn't do harm). If \code{recode=FALSE}, the
    categorical variables are assumed to be coded 1,2,3,...}
  \item{allout}{logical. If \code{TRUE}, the regular
    \code{flexmix}-output is given out for every single number of
    clusters, which can create a huge output object.}
  \item{control}{list of control parameters for \code{flexmix}, for
    details see the help page of \code{\link[flexmix]{FLXcontrol-class}}.}
  \item{silent}{logical. This is passed on to the
    \code{\link{try}}-function. If \code{FALSE}, error messages from
    failed runs of \code{flexmix} are suppressed. (The information that
    a \code{flexmix}-error occurred is still given out if
    \code{verbose=TRUE}).}
}
    


\details{
  Sometimes flexmix produces errors because of degenerating covariance
  matrices, too small clusters etc. \code{flexmixedruns} tolerates these
  and treats them as non-optimal runs. (Higher \code{simruns} or
  different \code{control} may be required to get a valid solution.)
  
  General documentation on flexmix can be found in 
  Friedrich Leisch's "FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R",
  \url{http://cran.r-project.org/web/packages/flexmix/vignettes/flexmix-intro.pdf}
}
 
\value{
  A list with components
  \item{optsummary}{summary object for \code{flexmix} object with
    optimal number of components.}
  \item{optimalk}{optimal number of components.}
  \item{errcount}{vector with numbers of EM runs for each number of
    components that led to flexmix errors.}
  \item{flexout}{if \code{allout=TRUE}, list of flexmix output objects
    for all numbers of components, for details see the help page of
    \code{\link[flexmix]{flexmix-class}}. Slots that can be used
    include for example \code{cluster} and \code{components}. So
    if \code{fo} is the \code{flexmixedruns}-output object,
    \code{fo$flexout[[fo$optimalk]]@cluster} gives a component number
    vector for the observations (maximum posterior rule), and
    \code{fo$flexout[[fo$optimalk]]@components} gives the estimated
    model parameters, which for \code{lcmixed} and therefore
    \code{flexmixedruns} are called
    \describe{
    \item{center}{mean vector}
    \item{cov}{covariance matrix}
    \item{pp}{list of categorical variable-wise category probabilities}
    }
    If \code{allout=FALSE}, only the flexmix output object for the
    optimal number of components, i.e., the \code{[[fo$optimalk]]}
    indexing above can then be omitted.
  }
  \item{bicvals}{vector of values of the BIC for each number of
    components.}
  \item{ppdim}{vector of categorical variable-wise numbers of
    categories.}
  \item{discretelevels}{list of levels of the categorical variables
    belonging to what is treated by \code{flexmixedruns} as category
    1, 2, 3 etc.} 
}

\references{
  Hennig, C. and Liao, T. (2010) Comparing latent class and
  dissimilarity based clustering for mixed type variables with
  application to social stratification. Research report no. 308,
  Department of Statistical Science, UCL.
  \url{http://www.ucl.ac.uk/Stats/research/reports/psfiles/rr308.pdf}
    Revised version accepted for publication
  by Journal of the Royal Statistical Society
  Series C.

}


\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{lcmixed}}, \code{\link[flexmix]{flexmix}},
  \code{\link[flexmix]{FLXcontrol-class}},
  \code{\link[flexmix]{flexmix-class}},
  \code{\link{discrete.recode}}.}


\examples{
  set.seed(776655)
  v1 <- rnorm(100)
  v2 <- rnorm(100)
  d1 <- sample(1:5,100,replace=TRUE)
  d2 <- sample(1:4,100,replace=TRUE)
  ldata <- cbind(v1,v2,d1,d2)
  fr <- flexmixedruns(ldata,
    continuous=2,discrete=2,simruns=2,n.cluster=2:3,allout=FALSE)
  print(fr$optimalk)
  print(fr$optsummary)
  print(fr$flexout@cluster)
  print(fr$flexout@components)
}

\keyword{cluster}% __ONLY ONE__ keyword per line
