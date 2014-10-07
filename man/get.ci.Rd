\name{get.ci}
\alias{get.ci}
\title{Compute the confidence interval of cellularity and ploidy}

\description{
  Compute the confidence intervals for a cellularity-ploidy posterior probability matrix.
}

\usage{
   get.ci(cp.table, level = 0.95)
}
\arguments{
  \item{cp.table}{list obtained with \code{\link{baf.model.fit}} or \code{\link{mufreq.model.fit}}.}
  \item{level}{Confidence level.}
}

\value{
  A list with 6 items:
  \item{values.ploidy}{matrix of ploidy values with respective posterior probability.}
  \item{confint.ploidy}{boundaries of the confidence interval of the estimated ploidy.}
  \item{max.ploidy}{point estimate of the ploidy value that has the maximum posterior probability.}
  \item{values.cellularity}{matrix of cellularity values with respective posterior probability.}
  \item{confint.cellularity}{boundaries of the confidence interval of the estimated cellularity.}
  \item{max.cellularity}{point estimate of the cellularity value that has the maximum posterior probability.}
}

\seealso{
  \code{\link{cp.plot}}.
}

\examples{
  data(CP.example)
  CP.example.ci <- get.ci(CP.example)
  str(CP.example.ci)
}
