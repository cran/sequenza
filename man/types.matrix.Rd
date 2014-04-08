\name{types.matrix}
\alias{types.matrix}
\title{Creates a dataframe of type tags}

\description{
  Type tags are a utensil to distinguish genomic positions by their copy number state and the number of mutated alleles. This function creates a dataframe of all possible type tags, given the copy number of the normal sample and the range of possible copy numbers in the tumor sample.
}

\usage{
  types.matrix(CNt.min, CNt.max, CNn = 2)
}

\arguments{
  \item{CNt.min}{minimum copy number in the tumor.}
  \item{CNt.max}{maximum copy number in the tumor.}
  \item{CNn}{copy number of the normal sample.}
}

\value{
\code{types.matrix} returns a dataframe with the following columns:
\item{CNn}{number of alleles in the normal sample.}
\item{CNt}{numbers of alleles in the aberrant sample}
\item{Mt}{number of mutated alleles}
}
\details{
    A type tag consists of 3 integers signifying the copy number in the normal and tumor samples and the number of mutated alleles. \code{types.matrix} returns all the possible type tags that can be created within the given range of tumor copy numbers.
}

\seealso{
  \code{theoretical.mufreq}, \code{theoretical.depth.ratio}, \code{theoretical.baf}, \code{model.points}.
}

\examples{
  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome diploid.
  types.matrix(CNt.min = 0, CNt.max = 4, CNn = 2 )

  ## Generate matrix types from 0 to 4 copy number, being the non-tumor chromosome monoploid.
  types.matrix(CNt.min = 0, CNt.max = 4, CNn = 1 )

}
