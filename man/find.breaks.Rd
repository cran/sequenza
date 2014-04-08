\name{find.breaks}
\alias{find.breaks}
\alias{segment.breaks}
\title{Segmentation of sequencing data using an allele-specific copy number algorithm}

\description{
  This function uses \code{\link{aspcf}} from the package \pkg{copynumber} to segment depth ratio and B-allele frequency obtained from sequencing data.
}

\usage{
find.breaks(seqz.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...)
segment.breaks(seqz.tab, breaks, min.reads.baf = 1, weighted.mean = TRUE)
}

\arguments{
  \item{seqz.baf}{an seqz file containing only the heterozygous positions.}
  \item{seqz.tab}{a complete seqz file.}
  \item{gamma, kmin, baf.thres, verbose}{arguments passed to \code{aspcf}.}
  \item{breaks}{breaks as output by \code{find.breaks}.}
  \item{min.reads.baf}{threshold on the depth of the positions included to calculate the average BAF for segment.}
  \item{weighted.mean}{boolean to select if the segments have to calculated using the read depth as a weights to calculate depth ratio and B-allele frequency means.}
  \item{...}{additional arguments passed to \code{\link{aspcf}}.}
}

\details{
  \pkg{copynumber} is a package to perform efficient segmentation of SNP-array data. The function \code{find.breaks} uses an allele-specific algorithm from the \pkg{copynumber} package to find break points, where the default parameters have been optimized for sequencing data, but a careful choice of an optimal \code{gamma} value is advised.
}

\examples{

  \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")
# read all the chromosomes:
seqz.data  <- read.seqz(data.file)
# Gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
seqz.data  <- read.seqz(data.file, chr.name = 12)

# Correct the coverage of the loaded chromosome:
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                           gc.vect[as.character(seqz.data$GC.percent)]
# Select the heterozygous positions
seqz.hom  <- seqz.data$zygosity.normal == 'hom'
seqz.het  <- seqz.data[!seqz.hom, ]
# Detect breakpoints
breaks <- find.breaks(seqz.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
segment.breaks(seqz.data, breaks = breaks)
   }
}
