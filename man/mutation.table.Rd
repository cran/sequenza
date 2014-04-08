\name{mutation.table}
\alias{mutation.table}
\title{Identify mutations}

\description{
 This function extracts positions from an seqz file that differ from the normal genome, applying various filters.
}

\usage{
mutation.table(seqz.tab, mufreq.treshold = 0.15, min.reads = 40, min.reads.normal = 10,
               max.mut.types = 3, min.type.freq = 0.9, min.fw.freq = 0, segments = NULL)
}

\arguments{
  \item{seqz.tab}{an seqz table, as output from \code{\link{read.seqz}}.}
  \item{mufreq.treshold}{mutation frequency threshold.}
  \item{min.reads}{minimum number of reads above the quality threshold to accept the mutation call.}
  \item{min.reads.normal}{minimum number of reads used to determine the genotype in the normal sample.}
  \item{max.mut.types}{maximum number of different base substitutions per position. Integer from 1 to 3 (since there are only 4 different bases). Default is 3, to accept \dQuote{noisy} mutation calls.}
  \item{min.type.freq}{minimum frequency of aberrant types.}
  \item{min.fw.freq}{minimum frequency of variant reads detected in the forward strand. Setting it to 0, all the variant calls with strand frequency in the interval outside 0 and 1, margin not comprised, would be discarded.}
  \item{segments}{if specified, the values of depth ratio would be taken from the segments rather than from the raw data.}
}

\details{
  Calling mutations in impure tumor samples is a difficult task, because the degree of contamination by normal cells affects the measured mutation frequency. In highly impure samples, where the normal cells comprise the major component of the sample, mutations can be so diluted that it can be difficult to distinguish them from sequencing errors.

  The function \code{mutation.table} tries to separate true mutations from sequencing errors, based on the given threshold. In samples with low contamination, it should even be possible to catch sub-clonal mutations using this function.
  
  This function identifes only those mutations occuring in positions that are homozygous in the normal genome.
}

\value{
A data frame, which in addition to some of the columns of the seqz table, contains the following two columns:
\item{F}{the mutation frequency}
\item{mutation}{a character representation of the mutation. For example, a mutation from \samp{A} in the normal to \samp{G} in the tumor is annotated as \samp{A>G}.}

}
\examples{

   \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")
seqz.data  <- read.seqz(data.file)

## Normalize coverage by GC-content
gc.stats <- gc.norm(x = seqz.data$depth.ratio,
                    gc = seqz.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                           gc.vect[as.character(seqz.data$GC.percent)]

## Extract mutations
mut.tab   <- mutation.table(seqz.data, mufreq.treshold = 0.15,
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9)
mut.tab <- na.exclude(mut.tab)
   }
}
