\name{model.points}

\alias{model.points}
\alias{baf.model.points}
\alias{mufreq.model.points}

\title{Generate B-allele frequency, mutation frequency and depth ratios
    at given model points, cellularity and ploidy values}

\description{
    The \code{baf.model.points} and \code{mufreq.model.points} functions
    combine \code{theoretical_baf}, \code{theoretical_mufreq} and
    \code{theoretical_depth_ratio} to model the theoretical respective values
        at known values of cellularity and ploidy.
}

\usage{
    baf.model.points(cellularity, ploidy, baf_types, avg.depth.ratio)
    mufreq.model.points(cellularity, ploidy, mufreq_types, avg.depth.ratio)
}

\arguments{
    \item{cellularity}{fraction of tumor cells in the sample.}
    \item{ploidy}{2 * ratio between total DNA content in a tumor cell
        and a normal cell.}
    \item{baf_types}{matrix with the sets of copy numbers and number of
        mutated alleles over which to model mutation frequency and depth ratio.
        The matrix can be generated with \code{\link{baf.types.matrix}}.}
    \item{mufreq_types}{matrix with the sets of copy numbers and number of
        mutated alleles over which to model mutation frequency and depth ratio.
        The matrix can be generated with \code{\link{mufreq.types.matrix}}.}
    \item{avg.depth.ratio}{average normalized depth ratio.}
}
\value{
    For \code{baf.model.points} a data.frame with two columns:
        \item{BAF}{modelled values of B-allele frequency.}
        \item{depth_ratio}{modelled values of depth ratio.}
    For \code{mufreq.model.points} a data.frame with two columns:
        \item{mufreqs}{modelled values of mutation frequency.}
        \item{depth_ratio}{modelled values of depth ratio.}
}

\details{
    The \code{baf.model.points} and \code{mufreq.model.points} functions
    generate the theoretical values of B-allele frequency, mutation frequency
    and depth ratio for the given type tags. To learn more about type tags
    see \code{\link{types.matrix}}.
}

\seealso{
    \code{\link{types.matrix}}, \code{\link{theoretical.depth.ratio}},
    \code{\link{theoretical.baf}} \code{\link{theoretical.mufreq}}.
}

\examples{
    # Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
    # and max = 4 on an originally diploid genome:
    types <- baf.types.matrix(CNt.min = 0, CNt.max = 4, CNn = 2)
    cbind(types, baf.model.points(cellularity = 0.5, ploidy = 2,
        baf_types = types, avg.depth.ratio = 1))
    # Simulate a cellularity of 0.5, ploidy of 2 and types from min CNt 0
    # and max = 4 on an originally monoallelic genome:
    types <- mufreq.types.matrix(CNt.min = 0, CNt.max = 4, CNn = 1)
    cbind(types, mufreq.model.points(cellularity = 0.5, ploidy = 2,
        mufreq_types = types, avg.depth.ratio = 1))
}
