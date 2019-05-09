\name{types.matrix}
\alias{types.matrix}
\alias{baf.types.matrix}
\alias{mufreq.types.matrix}
\title{Creates a matrix of type tags}

\description{
    Type tags are a utensil to distinguish genomic positions by their copy
    number state, number A and B alleles and the number of mutated alleles.
    This function creates a matrix of all possible type tags, given the copy
    number of the normal sample and the range of possible copy numbers in
    the tumor sample.
}

\usage{
    baf.types.matrix(CNt.min, CNt.max, CNn = 2)
    mufreq.types.matrix(CNt.min, CNt.max, CNn = 2)
}

\arguments{
    \item{CNt.min}{minimum copy number in the tumor.}
    \item{CNt.max}{maximum copy number in the tumor.}
    \item{CNn}{copy number of the normal sample.}
}

\value{
    \code{baf.types.matrix} returns a data.frame with the 3 columns:
        \item{CNn}{number of alleles in the normal sample.}
        \item{CNt}{numbers of alleles in the tumor sample.}
        \item{B}{number of B alleles in the tumor sample.}

    \code{mufreq.types.matrix} returns a data.frame with the 3 columns:
        \item{CNn}{number of alleles in the normal sample.}
        \item{CNt}{numbers of alleles in the tumor sample.}
        \item{Mt}{number of mutated alleles in the tumor sample.}
}

\details{
    A type consists of 3 integers signifying the copy number in the normal
    and tumor samples and the number of B alleles (\code{baf.types.matrix}) or
    mutated alleles (\code{mufreq.types.matrix}).
    The two functions return all the possible types combination within the
    range of tumor copy numbers in the arguments (\code{CNt.min:CNt.max}).
}

\seealso{
    \code{theoretical_mufreq}, \code{theoretical_depth_ratio},
    \code{theoretical_baf}, \code{model_points}.
}

\examples{
    ## Generate matrix types from 0 to 4 copy number, being the
    ## non-tumor chromosome diploid.
    baf.types.matrix(CNt.min = 0, CNt.max = 4, CNn = 2 )

    ## Generate matrix types from 0 to 4 copy number, being the
    ## non-tumor chromosome monoploid.
    mufreq.types.matrix(CNt.min = 0, CNt.max = 4, CNn = 1 )
}
