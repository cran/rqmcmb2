\name{rqmcmb.plot}
\alias{rqmcmb.plot}
\title{Plot the MCMB Sequences}

\description{
  Plots the MCMB sequences. The dotted line is the initial estimate of the
  parameter vector from \code{rq}.  The two solid lines are the CI limits of
  size alpha centered at the MCMB parameter estimates}


\usage{
  rqmcmb.plot(obj, alpha=.10)
}

\arguments{
  \item{obj}{The output of rqmcmb() function}
  \item{alpha}{Confidence level}
}

\keyword{regression}
\keyword{hplot}
