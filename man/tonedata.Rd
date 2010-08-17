\name{tonedata}
\alias{tonedata}
\docType{data}
\title{Tone perception data}
\description{
  The tone perception data stem
  from an experiment of Cohen (1980) and have been analyzed in de Veaux
  (1989).
  A pure fundamental tone was played to a
  trained musician. Electronically generated overtones were added, determined 
  by a stretching ratio of \code{stretchratio}. \code{stretchratio=2.0}
  corresponds to the harmonic pattern
  usually heard in traditional definite pitched instruments. The musician was
  asked to tune an adjustable tone to the octave above the fundamental tone.
  \code{tuned} gives the ratio of the adjusted tone to the fundamental,
  i.e. \code{tuned=2.0} would be the correct tuning for all
  \code{stretchratio}-values.
  The data analyzed here belong to 150 trials 
  with the same musician. In the original study, there were four further
  musicians. 
}
\usage{data(tonedata)}
\format{A data frame with 2 variables \code{stretchratio} and
  \code{tuned} and 150 cases.}
\source{Cohen, E. A. (1980) \emph{Inharmonic tone
    perception}. Unpublished Ph.D. dissertation, Stanford University}
\references{
  de Veaux, R. D. (1989) Mixtures of Linear Regressions,
  \emph{Computational Statistics and Data Analysis} 8, 227-245.
}
\keyword{datasets}
