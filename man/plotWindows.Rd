\name{plotWindows}
\alias{plotWindows}
\title{Plot a windowed chromosome}

\description{
  \code{plotWindows} is a convenient plot function that visualizes a specified chromosome attribute using the data frame produced by the function \code{windowValues}.
}

\usage{
plotWindows(seqz.window, m.lty = 1, m.lwd = 3,
             m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
             n.min = 1, xlim, ylim, add = FALSE, ...)
}

\arguments{
  \item{seqz.window}{data frame of base-pair windows and corresponding statistics to be plotted. A list of such data frames can be output from \code{\link{windowValues}}}
  \item{m.lty}{line type used for plotting mean values.}
  \item{m.lwd}{line width used for plotting mean values.}
  \item{m.col}{line color used for plotting mean values.}
  \item{q.bg}{background color for the area between the 0.25 and 0.75 quartiles.}
  \item{log2.plot}{logical, if TRUE values are log2 scaled.}
  \item{n.min}{minimum number of data points required for a binned window to be plotted.}
  \item{xlim}{limits of the x axis.}
  \item{ylim}{limits of the y axis.}
  \item{add}{logical, if TRUE the plot will be added to an existing opened device.}
  \item{...}{any other arguments accepted by \code{\link{plot}}.}
}

\seealso{
  \code{\link{chromosome.view}},
}

\examples{

   \dontrun{

data.file <-  system.file("data", "example.seqz.txt.gz", package = "sequenza")
seqz.data <- read.seqz(data.file)
# 1Mb windows, each window is overlapping with 1 other adjacent window: depth ratio
seqz.ratio <- windowValues(x = seqz.data$depth.ratio,
                           positions = seqz.data$position,
                           chromosomes = seqz.data$chromosome, window = 1e6,
                           weight = seqz.data$depth.normal, start.coord = 1, overlap = 1)

plotWindows(seqz.ratio[[1]], log2.plot = FALSE,
            ylab = "Depth ratio", xlab = "Position (bases)",
            main = names(seqz.ratio)[1], las = 1, n.min = 1,
            ylim = c(0, 2.5))

plotWindows(seqz.ratio[[17]], log2.plot = FALSE,
            ylab = "Depth ratio", xlab = "Position (bases)",
            main = names(seqz.ratio)[1], las = 1, n.min = 1,
            ylim = c(0, 2.5))
   }
}
