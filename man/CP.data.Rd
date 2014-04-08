\name{CP.example}
\docType{data}
\alias{CPsnp.example}
\alias{CP.example}

\title{Example of cellularity and ploidy results}

\description{
  Examples of objects resulting from the Bayesian estimation over a set of cellularity and ploidy parameters, as returned by \code{\link{baf.model.fit}} or \code{\link{mufreq.model.fit}}.
}

\usage{
  data(CP.example)
  data(CPsnp.example)
}

\format{
  A list containing three items:
 \describe{ 
   \item{ploidy}{numeric vector of tested ploidy values.}
   \item{cellularity}{numeric vector of tested cellularity values.}
   \item{loglik}{numeric matrix of log-likelihoods for each \var{(ploidy, cellularity)} pair.}
  }
}

\examples{
  data(CP.example)
  str(CP.example)

  data(CPsnp.example)
  str(CPsnp.example)

  ## Visualization of the object
  image(x = CP.example$ploidy,
        y = CP.example$cellularity,
        z = CP.example$loglik)

  ## A better plot
  cp.plot(CP.example)
  cp.plot.contours(CP.example, add = TRUE)
}

\keyword{datasets}
