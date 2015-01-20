\name{VarScan2seqz}
\alias{VarScan2seqz}
\title{Import VarScan output files, for use with the Sequenza algorithm.}

\description{
   This function reads the output files from the \command{somatic} and \command{copynumber} programs, which are part of VarScan version 2. The resulting object can be used, in lieu of a \dQuote{seqz} object, within Sequenza.
}

\usage{
  VarScan2seqz(varscan.somatic, varscan.copynumber = NULL, normal_var_freq = 0.25)
}

\arguments{
  \item{varscan.somatic}{a data frame from the output file of \command{somatic}.}
  \item{varscan.copynumber}{a data frame from the output file of \command{copynumber} (optional).}
  \item{normal_var_freq}{defines the variant allele frequency threshold on the normal genotype to filter out false positive heterozygous position.}
}

\details{
  It is possible to run sequenza with the output of VarScan2 version 2.3.6. Earlier versions may also work but have not been tested.
}

\value{
A data frame with the seqz format, as returned by \code{\link{read.seqz}}.
}

\references{
\url{http://varscan.sourceforge.net}
}
