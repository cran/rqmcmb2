\name{scoref}
\alias{scoref}
\title{ Score Function for Quantile Regression}
\description{
  The score function returns tau if the residual r is positive, tau-1 if 
  the residual r is negative, and 0 if the residual r is equal to 0}

\usage{
  scoref(r, tau)
}

\arguments{
  \item{r}{The residual}
  \item{tau}{The desire quantile}
}

\references{
   Koenker, R. W. and Bassett, G. W. (1978) Regression
     quantiles.  Econometrica, 46, 33-50.}

\keyword{internal}


