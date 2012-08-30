\name{discrete.recode}
\alias{discrete.recode}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Recodes mixed variables dataset}
\description{
  Recodes a dataset with mixed continuous and categorical variables so
  that the continuous variables come first and the categorical variables
  have standard coding 1, 2, 3,... (in lexicographical ordering of
  values coerced to strings).
}
\usage{
  discrete.recode(x,xvarsorted=TRUE,continuous=0,discrete)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame. The data need to be organised
    case-wise, i.e., if there are categorical variables only, and 15
    cases with values c(1,1,2) on the 3 variables, the data matrix needs
    15 rows with values 1 1 2. (Categorical variables could take numbers
    or strings or anything that can be coerced to factor levels as values.)}
  \item{xvarsorted}{logical. If \code{TRUE}, the continuous variables
    are assumed to be the first ones, and the categorical variables to
    be behind them.}
  \item{continuous}{vector of integers giving positions of the
    continuous variables. If \code{xvarsorted=TRUE}, a single integer,
    number of continuous variables.}
  \item{discrete}{vector of integers giving positions of the
    categorical variables (the variables need to be coded in such a way that
    \code{\link{data.matrix}} converts them to something numeric). If
    \code{xvarsorted=TRUE}, a single integer, number of categorical variables.}
}

\value{
  A list with components
  \item{data}{data matrix with continuous variables first and
  categorical variables in standard coding behind them.}
  \item{ppdim}{vector of categorical variable-wise numbers of
    categories.}
  \item{discretelevels}{list of levels of the categorical variables
    belonging to what is treated by \code{flexmixedruns} as category
    1, 2, 3 etc.} 
  \item{continuous}{number of continuous variables.}
  \item{discrete}{number of categorical variables.}
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{lcmixed}}}

\examples{
  set.seed(776655)
  v1 <- rnorm(20)
  v2 <- rnorm(20)
  d1 <- sample(c(2,4,6,8),20,replace=TRUE)
  d2 <- sample(1:4,20,replace=TRUE)
  ldata <- cbind(v1,d1,v2,d2)
  lc <-
  discrete.recode(ldata,xvarsorted=FALSE,continuous=c(1,3),discrete=c(2,4))
  data(Cars93)
  Cars934 <- Cars93[,c(3,5,8,10)]
  cc <- discrete.recode(Cars934,xvarsorted=FALSE,continuous=c(2,3),discrete=c(1,4))
}

\keyword{manip}% __ONLY ONE__ keyword per line
