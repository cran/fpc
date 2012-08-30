\name{cat2bin}
\alias{cat2bin}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Recode nominal variables to binary variables}
\description{
  Recodes a dataset with nominal variables so that the nominal
  variables are replaced by binary variables for the categories.
}
\usage{
  cat2bin(x,categorical=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{data matrix or data frame. The data need to be organised
    case-wise, i.e., if there are categorical variables only, and 15
    cases with values c(1,1,2) on the 3 variables, the data matrix needs
    15 rows with values 1 1 2. (Categorical variables could take numbers
    or strings or anything that can be coerced to factor levels as values.)}
  \item{categorical}{vector of numbers of variables to be recoded.}
}

\value{
  A list with components
  \item{data}{data matrix with variables specified in \code{categorical}
    replaced by 0-1 variables, one for each category.}
  \item{variableinfo}{list of lists. One list for every variable in the
    original dataset, with four components each, namely \code{type}
    (\code{"categorical"} or \code{"not recoded"}), \code{levels}
  (levels of nominal recoded variables in order of binary variable in
  output dataset), \code{ncat} (number of categories for recoded
  variables), \code{varnum} (number of variables in output dataset
  belonging to this original variable).}  
}

\author{Christian Hennig
  \email{c.hennig@ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{discrete.recode}}}

\examples{
  set.seed(776655)
  v1 <- rnorm(20)
  v2 <- rnorm(20)
  d1 <- sample(1:5,20,replace=TRUE)
  d2 <- sample(1:4,20,replace=TRUE)
  ldata <- cbind(v1,v2,d1,d2)
  lc <- cat2bin(ldata,categorical=3:4)
}



\keyword{manip}% __ONLY ONE__ keyword per line
