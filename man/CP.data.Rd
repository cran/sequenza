\name{CP.example}
\alias{CPsnp.example}
\alias{CP.example}
\alias{CP}
\alias{CP.snp}

\title{Example of cellularity and ploidy results.}

\description{
  Examples of objects resulting from the Bayesian estimation over a set of cellularity and ploidy parameters.
}

\value{
 \item{x}{vector containing the tested ploidy values, sorted from the smaller to the bigger value.}
 \item{y}{vector containing the tested cellularity values, sorted from the smaller to the bigger value.}
 \item{z}{matrix with x columns and y rows, containing the log-likelihood of each tested ploidy and cellularity cobination.}
}

\examples{
  data(CP.example)
  str(CP.example)
  data(CPsnp.example)
  str(CPsnp.example)
  # Visualization of the object
  image(CP.example)
  # A better plot
  cp.plot(CP.example)
  cp.plot.contours(CP.example, add = TRUE)
}
