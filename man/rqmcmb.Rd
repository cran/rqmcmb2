\name{rqmcmb}
\alias{rqmcmb}


\title{Markov Chain Marginal Bootstrap for Quantile
Regression} 


\description{MCMB for Quantile Regression (also see
quantreg package by Roger Koenker)}



\usage{
rqmcmb(x=x, y=y, tau=0.5, K=100, int=TRUE,
plotTheta=FALSE)
}


\arguments{
  \item{x}{a data matrix (n by p) for the design variables whose rows correspond to cases}
  \item{y}{a response vector of length n}
  \item{tau}{a percentile level between 0 and 1. Default at 0.5 for the median}
  \item{K}{length of the MCMB sequence. Default is 100}
  \item{plotTheta}{TRUE or FALSE for plotting the MCMB sequence. Default to FALSE} 
  \item{int}{should be set to TRUE if the intercept is to be included in the model, and 
to FALSE if no intercept is desired. Default is TRUE.} 
}


\value{
  A list with the following components:
  \item{coef}{the parameter estimate from rq()}
  \item{theta}{a matrix containing the MCMB sequence. The first row is  the initial 
parameter estimate from rq()}
  \item{success}{returns 1 if MCMB is successful.  A value of 0 
indicates that the program fails to return a desired MCMB sequence}
  \item{cn}{condition number of the X'X matrix.}
}



\references{
Kocherginsky, M. (2003). Extensions of Markov
Chain Marginal Bootstrap. Ph.D Thesis, University of Illinois
Urbana-Champaign.
 
Kocherginsky, M., He, X. and Mu, Y. (2003).
Practical Confidence Intervals for Regression Quantiles. Preprint.

He, X. and Hu, F. (2002). Markov Chain Marginal Bootstrap. 
Journal of the American Statistical Association , Vol. 97, no. 459, 
783-795. 
}


\author{ Maria Kocherginsky
(mkocherg@health.bsd.uchicago.edu) and Xuming He (x-he@uiuc.edu)}


\section{WARNING }{The MCMB may not be suitable for
problems of small sample sizes. Severe collinearity in the x matrix
could also be harmful.} 



\seealso{\code{\link[quantreg]{rq}}}


\examples{
  library(quantreg)
  
  x <- cbind(rnorm(100), runif(100))
  y <- rnorm(100)

  #generate the MCMB sequence:
  mcmb <- rqmcmb(x, y, tau=.5, plotTheta=FALSE)

  #get MCMB estimates of mean, SD, and CI:
  rqmcmb.ci(mcmb)

  #plot the MCMB sequences:
  rqmcmb.plot(mcmb)
  

}

\keyword{regression}
\keyword{robust}
