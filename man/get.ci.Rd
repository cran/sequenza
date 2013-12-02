\name{get.ci}
\alias{get.ci}
\title{Compute the confidence interval of cellularity and ploidy}

\description{
  Compute the confidence intervals for a cellularity-ploidy likelihood matrix.
}

\usage{
   get.ci(cp.table, interval = 0.95)
}
\arguments{
  \item{cp.table}{list obtained with \code{baf.model.fit} or \code{mufreq.model.fit}.}
  \item{interval}{size of confidence interval (in decimal). Defaults to 0.95}
}

\value{
  \item{values.x}{matrix of ploidy values with respective likelihood.}
  \item{confint.x}{boundaries of the confidence interval of the estimated ploidy.}
  \item{max.x}{point estimate of the ploidy value that has the maximum likelihood.}
  \item{values.y}{matrix of cellularity values with respective likelihood.}
  \item{confint.y}{boundaries of the confidence interval of the estimated cellularity.}
  \item{max.y}{point estimate of the cellularity value that has the maximum likelihood.}
}

\seealso{
  \code{baf.model.fit}, \code{mufreq.model.fit}.
}