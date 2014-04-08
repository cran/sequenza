\name{model.points}
\alias{model.points}
\title{Generate mutation frequency and depth ratio model point at given cellularity and ploidy value}

\description{
  \code{model.points} combines \code{theoretical.mufreq} and \code{theoretical.depth.ratio} to model the theoretical mutation frequency and depth ratio based on given values of cellularity and ploidy.
}

\usage{
model.points(cellularity, ploidy, types, avg.depth.ratio)
}

\arguments{
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{2 * ratio between total DNA content in a tumor cell and a normal cell.}
  \item{types}{matrix with the sets of copy numbers and number of mutated alleles over which to model mutation frequency and depth ratio. The matrix can be generated with \code{\link{types.matrix}}.}
  \item{avg.depth.ratio}{average normalized depth ratio.}
}
\value{
A matrix of two columns:
  \item{mufreqs}{modelled values of mutation frequency for the given types.}
  \item{depth.ratio}{modelled values of depth ratio for the given types.}
}

\details{
  \code{model.points} generates the theoretical values of mutation frequency and depth ratio for the given type tags. To learn more about type tags see \code{\link{types.matrix}}.
}

\seealso{
  \code{\link{types.matrix}}, \code{\link{theoretical.depth.ratio}}, \code{\link{theoretical.mufreq}}.
}

\examples{

# Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
# and max = 4 on an originally diploid genome:
types <- types.matrix(CNt.min = 0, CNt.max = 4, CNn = 2)
cbind(types, model.points(cellularity = 0.5, ploidy = 2,
                          types = types, avg.depth.ratio = 1))
# Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
# and max = 4 on an originally monoallelic genome:
types <- types.matrix(CNt.min = 0, CNt.max = 4, CNn = 1)
cbind(types, model.points(cellularity = 0.5, ploidy = 2,
                          types = types, avg.depth.ratio = 1))

}
