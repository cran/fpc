\name{weightplots}
\alias{weightplots}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ordered posterior plots}
\description{
  Ordered posterior plots for Gaussian mixture components, see Hennig (2010).
}
\usage{
  weightplots(mclustsummary, clusternumbers="all", clustercol=2,
                        allcol=grey(0.2+((1:ncol(mclustsummary$z))-1)*
                          0.6/(ncol(mclustsummary$z)-1)),
                        lty=rep(1,ncol(mclustsummary$z)),clusterlwd=3,
                        legendposition="none",
                        weightcutoff=0.01,ask=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mclustsummary}{output object from
    \code{\link[mclust]{summary.mclustBIC}}.}
  \item{clusternumbers}{\code{"all"} or vector of integers. Numbers of
    components for which plots are drawn.}
  \item{clustercol}{colour used for the main components for which a
    plot is drawn.}
  \item{allcol}{colours used for respective other components in plots in
    which they are not main components.}
  \item{lty}{line types for components.}
  \item{clusterlwd}{numeric. Line width for main component.}
  \item{legendposition}{\code{"none"} or vector with two coordinates in
    the plot, where a legend should be printed.}
  \item{weightcutoff}{numeric between 0 and 1. Observations are only taken
    into account for which the posterior probability for the main
    component is larger than this.}
  \item{ask}{logical. If \code{TRUE}, it sets \code{par(ask=TRUE)} in
    the beginning and \code{par(ask=FALSE)} after all plots were showed.}
  \item{...}{further parameters to be passed on to \code{\link{legend}}.}
}

\value{
  Invisible matrix of posterior probabilities \code{z} from
  \code{mclustsummary}.
}

\details{
  Shows posterior probabilities for observations belonging to all
  mixture components on the y-axis, with points ordered by posterior
  probability for main component.
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
  data(crabs)
  dc <- crabs[,4:8]
  cm <- mclustBIC(crabs[,4:8],G=9,modelNames="EEE")
  scm <- summary(cm,crabs[,4:8])
  weightplots(scm,clusternumbers=1:3,ask=FALSE)
  weightplots(scm,clusternumbers=1:3,allcol=1:9, ask=FALSE,
              legendposition=c(5,0.7))
# Remove ask=FALSE to have time to watch the plots.
}
\keyword{multivariate}
\keyword{cluster}


