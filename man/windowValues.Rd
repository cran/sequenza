\name{windowValues}
\alias{windowValues}
\alias{windowBf}
\title{Bins sequencing data for plotting}

\description{
  Given a variable with corresponding genomic positions, the function bins the
  values in windows of a specified size and calculates weighted mean and 25th
  and 75th percentile for each window. The resulting object are visualized by
  the function \code{plotWindows}.
}

\usage{
    windowValues(x, positions, chromosomes, window = 1e6, overlap = 0,
        weight = rep.int( x = 1, times = length(x)), start.coord = 1)
    windowBf(Af, Bf, good.reads, positions, chromosomes, window = 1e6,
        overlap = 0, start.coord = 1, conf = 0.95)
}
\arguments{
    \item{x}{variable to be windowed.}
    \item{positions}{base-pair positions.}
    \item{chromosomes}{names or numbers of the chromosomes.}
    \item{window}{size of windows used for binning data. Smaller windows will
        take more time to compute.}
    \item{overlap}{integer defining the number of overlapping windows. Default
        is 0, no overlap.}
    \item{weight}{weights to be assigned to each value of \code{x}, usually
        related to the read depth.}
    \item{start.coord}{coordinate at which to start computing the windows. If
        NULL, will start at the first position available.}
    \item{Af}{A-allele frequency for the \code{Bf} calculation.}
    \item{Bf}{B-allele frequency for the \code{Bf} calculation.}
    \item{good.reads}{number of reads passing filter for the \code{Bf}
        calculation.}
    \item{conf}{confidence intervals of the binned \code{Bf} value.}
}

\details{
    DNA sequencing produces an amount of data too large to be handled by
    standard graphical devices. In addition, for samples analyzed with older
    machines and with low or middle coverage (20x to 50x), measures such as
    read depth are subject to big variations due to technical noise.
    Using \code{windowValues} prior to plotting reduces the noise and the
    amount of data to be plotted.

    The binning of the B-allele frequency requires a separate function,
    \code{windowBf}, as the B-allele frequency calculation uses multiple
    values: \code{Af}, \code{Bf} and \code{good.reads}.

    The output of \code{windowValues} and \code{windowBf} can be used as input
    for \code{\link{plotWindows}}.
}

\value{
    a list of data.frame, one per chromosome. Each data.frame contains
    base-pair windows covering the chromosome. Each row of the data.frame
    correspond to a window and its weighted mean, 25th and 75th percentiles of
    the input values, and the number of data points within each window.
}

\seealso{
    \code{plotWindows}
}

\examples{

\dontrun{
    data.file <-  system.file("extdata", "example.seqz.txt.gz",
        package = "sequenza")
    seqz.data <- read.seqz(data.file)
    # 1Mb windows, each window is overlapping with 1 other
    # adjacent window: depth ratio
    seqz.ratio <- windowValues(x = seqz.data$depth.ratio,
        positions = seqz.data$position, chromosomes = seqz.data$chromosome,
        window = 1e6, weight = seqz.data$depth.normal, start.coord = 1,
        overlap = 1)

    seqz.hom  <- seqz.data$zygosity.normal == 'hom'
    seqz.het  <- seqz.data[!seqz.hom, ]
    # 1Mb windows, each window is overlapping with 1 other adjacent window:
    # B-allele frequency
    seqz.bafs  <- windowValues(x = seqz.het$Bf, positions = seqz.het$position,
        chromosomes = seqz.het$chromosome, window = 1e6,
        weight = seqz.het$depth.tumor, start.coord = 1, overlap = 1)
    # Repeat the same operation using windowBf
    seqz.bafs  <- windowBf(Af = seqz.het$Bf, Bf = seqz.het$Bf,
        good.reads = seqz.het$good.reads, positions = seqz.het$position,
        chromosomes = seqz.het$chromosome, window = 1e6,
        start.coord = 1, overlap = 1, conf = 0.95)
}
}
