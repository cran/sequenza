\name{windowValues}
\alias{windowValues}
\title{Bins sequencing data for plotting}

\description{
  Given a variable with corresponding genomic positions, the function bins the values in windows of a specified size and calculates weighted mean and 25th and 75th percentile for each window, which is useful for plotting with \code{plotWindows}.
}

\usage{
windowValues(x, positions, chromosomes, window = 1e6, overlap = 0, verbose = TRUE,
             weight = rep.int( x = 1, times = length(x)), start.coord = 1)
}
\arguments{
  \item{x}{variable to be windowed.}
  \item{positions}{base-pair positions.}
  \item{chromosomes}{names or numbers of the chromosomes.}
  \item{window}{size of windows used for binning data. Smaller windows will take more time to compute.}
  \item{overlap}{integer defining the number of overlapping windows. Default is 0, no overlap.}
  \item{verbose}{logical, indicating whether to print information about the chromosome being processed.}
  \item{weight}{weights to be assigned to each value of \code{x}, usually related to the read depth.}
  \item{start.coord}{coordinate at which to start computing the windows. If NULL, will start at the first position available.}
}

\details{
  DNA sequencing produces an amount of data too large to be handled by standard graphical devices. In addition, for samples analyzed with older machines and with low or middle coverage (20x to 50x), measures such as read depth are subject to big variations due to technical noise. Using \code{windowValues} prior to plotting reduces the amount of data that will be plotted, and reduces noise.

  The output of \code{windowValues} can be used as input for \code{plotWindows}.
}

\value{
  a list of matrices, one per chromosome. Each matrix contains base-pair windows covering the chromosome, and for each window the weighted mean, 25th and 75th percentiles of the input values, and the number of data points within each window.
}

\seealso{
  \code{plotWindows}, \code{mutation.table}
}

\examples{
   \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
abf.data <- read.abfreq(data.file)
# 1Mb windows, each window is overlapping with 1 other adjacent window: depth ratio
abf.ratio <- windowValues(x = abf.data$depth.ratio,
                           positions = abf.data$n.base,
                           chromosomes = abf.data$chromosome, window = 1e6,
                           weight = abf.data$depth.normal, start.coord = 1, overlap = 1)
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# 1Mb windows, each window is overlapping with 1 other adjacent window: B-allele frequency
abf.bafs  <- windowValues(x = abf.het$Bf,
                           positions = abf.het$n.base,
                           chromosomes = abf.het$chromosome, window = 1e6,
                           weight = abf.het$depth.sample,
                           start.coord = 1, overlap = 1)
   }
}
