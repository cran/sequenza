\name{gc.sample.stats}
\alias{gc.sample.stats}
\alias{gc.summary.plot}
\alias{mean_gc}
\alias{median_gc}


\title{Collect display and correct GC-content related coverage bias}

\description{
  Collect information and perform statistics of depth of coverage in
  relation with GC-content.
}

\usage{
  gc.sample.stats(file, col_types = "c--dd----d----", buffer = 33554432,
                 parallel = 2L, verbose = TRUE)
  gc.summary.plot(gc_list, mean.col = 1, median.col = 2,
     scale.subset = 1.5, ...)
  mean_gc(gc_list)
  median_gc(gc_list)
}

\arguments{
  \item{file}{name of a file in the seqz format.}
  \item{col_types}{a string describing the classes of each columns of the
        input file (see \code{\link{read_tsv}}). The default value corresponds
        to the columns of a seqz file used for carculating GC statistics.}
  \item{buffer}{maximal size of each chunk in bytes(see
        \code{\link{chunk.apply}}).}
  \item{parallel}{integer, number of threads used to process a seqz file
        (see \code{\link{chunk.apply}}).}
  \item{verbose}{logical. If TRUE (the default) the function retuns information
        in the console.}
  \item{gc_list}{a normal or tumor list resulting from the
     \code{\link{gc.sample.stats}} function.}
  \item{mean.col}{color for the mean in the summary plot.}
  \item{median.col}{color for the median in the summary plot.}
  \item{scale.subset}{scale the depth values to sho in the plot. A value of 1
    will show the average depth at the center of the plot.}
  \item{...}{additional parametrers from \code{\link{colorgram}}.}
}

\details{
  \code{gc.sample.stats} extracts depths and GC-content inforation for the
  tumor and the control samples from an seqz file
  it returns a list with 3 elements: \code{file.metrics}, \code{normal} and
  \code{tumor}.

  \code{file.metrics} is a \code{data.frame} serving as index of the seqz
  file; the \code{normal} and \code{tumor} objects contains each 3 ojects:
  \code{gc}, \code{depth} and \code{n}.

  \code{gc} and \code{depth} are vectors containing the recorded values of,
  respectively, GC and coverage depth. the \code{n} object is a matrix
  \code{gc}x\code{depth}, recording the number of time a certain
  \code{gc}/\code{depth} pairs is observed in the data.

}

\value{
A list with the following elements:
  \item{file.metrics}{index of the seqz file.}
  \item{tumor}{GC and coverage depth observations in the tumor sample.}
  \item{normal}{GC and coverage depth observations in the control sample.}
}


\examples{

\dontrun{

data.file <-  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")
# read all the chromosomes:
gc_info <- gc.sample.stats(data.file)

# mean values of depth coverage vs GC content

mean_gc(gc_info$normal)

# plot the information for the tumor and normal samples
par(mfrow=c(1, 2))
gc.summary.plot(gc_info$normal, main = "Normal GC stats")
gc.summary.plot(gc_info$tumor, main = "Tumor GC stats")
}
}