\name{find.breaks}
\alias{find.breaks}
\alias{segment.breaks}
\title{Segmentation of sequencing data using an allele-specific copy number algorithm}

\description{
  This function uses \code{\link{aspcf}} from the package \code{copynumber} to segment depth ratio and B-allele frequency obtained from sequencing data.
}

\usage{
find.breaks(abf.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...)
segment.breaks(abf.tab, breaks) 
}

\arguments{
  \item{abf.baf}{an ABfreq file containing only the heterozygous positions.}
  \item{abf.tab}{a complete ABfreq file.}  
  \item{gamma, kmin, baf.thres, verbose}{arguments passed to \code{aspcf}.}
  \item{breaks}{breaks as output \code{find.breaks}.}
  \item{...}{additional arguments passed to \code{aspcf}.}
}

\details{
  \code{copynumber} is a package to perform efficient segmentation of SNP-array data. The function \code{find.breaks} uses an allele-specific algorithm from the \code{copynumber} package to find break points, where the default parameters have been optimized for sequencing data, but a careful choice of an optimal \code{gamma} value is advised.
}

\examples{

  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
# read all the chromosomes:
abf.data  <- read.abfreq(data.file)
# Gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
abf.data  <- read.abfreq(data.file, chr.name = 12)

# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]
# Select the heterozygous positions
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# Detect breakpoints
breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
segment.breaks(abf.data, breaks = breaks)
   }
}
