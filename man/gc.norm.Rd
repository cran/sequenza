\name{gc.sample.stats}
\alias{gc.sample.stats}
\alias{gc.norm}

\title{Normalize depth ratio values for GC-content bias}

\description{
  Detects and bias in the depth ratio values driven by varying GC-content.
}

\usage{
gc.sample.stats(file, gz = TRUE)
gc.norm(x, gc)
}

\arguments{
  \item{file}{name of a file in the seqz format.}
  \item{x}{vector of values to be normalized by GC-content, typically depth ratio values.}
  \item{gc}{vector of relative GC-content values for x.}
  \item{gz}{logical. If TRUE (the default) the function expects a gzipped file.}
}

\details{
  \code{gc.norm} detects bias in \code{x} driven by \code{gc}. Specifically, for each value of \code{gc}, summary statistics are calculated for the corresponding values of \code{x}. These statistics can then be used to normalize \code{x} for \code{gc}.

  \code{gc.sample.stats} extracts depth ratio and GC-content from an seqz file, and then uses \code{gc.norm} on the results.

}

\value{
A list with the following elements:
  \item{raw}{quartiles of \code{x} for each value of \code{gc}}
  \item{adj}{median-normalized values of \code{raw}}
  \item{gc.values}{vector of different GC-content values observed}
  \item{raw.mean}{mean of \code{x} for each value of \code{gc}}
  \item{raw.median}{median \code{x} for each value of \code{gc}}
  \item{file.metrics}{only from \code{gc.sample.stats}.}
}


\examples{

  \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")
# read all the chromosomes:
seqz.data  <- read.seqz(data.file)
# Normalize coverage by GC-content
gc.stats <- gc.norm(x = seqz.data$depth.ratio,
                    gc = seqz.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                           gc.vect[as.character(seqz.data$GC.percent)]

# Alternatively gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
seqz.data  <- read.seqz(data.file, chr.name = 12)
# Correct the coverage of the loaded chromosome:
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                           gc.vect[as.character(seqz.data$GC.percent)]

   }
}
