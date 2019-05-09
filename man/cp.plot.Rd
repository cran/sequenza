\name{cp.plot}
\alias{cp.plot}
\alias{cp.plot.contours}
\alias{get.ci}

\title{Plot log-posterior probability for the output of the
 \code{\link{sequenza.fit}} function}

\description{
  This function uses the \code{\link{colorgram}} function from the package
  \pkg{squash} to plot log-posterior probability for the tested combinations of
  cellularity and ploidy
}

\usage{
   cp.plot(cp.table, xlab = "Ploidy", ylab = "Cellularity",
           zlab = "Scaled rank LPP",
           colFn = colorRampPalette(c('white', 'lightblue')), ...)
   cp.plot.contours(cp.table, likThresh = c(0.95), alternative = TRUE,
                    col = palette(), legend.pos = "bottomright", pch = 18,
                    alt.pch = 3, ...)
   get.ci(cp.table, level = 0.95)
}

\arguments{
  \item{cp.table}{list, as output from \code{\link{baf.model.fit}} or \code{\link{mufreq.model.fit}}.}
  \item{xlab}{xlab parameter as in the function \code{\link{colorgram}}.}
  \item{ylab}{ylab parameter as in the function \code{\link{colorgram}}.}
  \item{zlab}{zlab parameter as in the function \code{\link{colorgram}}.}
  \item{colFn}{colFn parameter as in the function \code{\link{colorgram}}.}
  \item{likThresh}{vector of quantiles to define tresholds for the confindent regions.}
  \item{alternative}{boolean parameter, if \code{TRUE} the alternative solutions are computed and plotted.}
  \item{col}{vector of colors.}
  \item{legend.pos}{position for placing the legend.}
  \item{pch}{character used to indicate the point estimate.}
  \item{alt.pch}{if \code{alternative} is set to \code{TRUE} defines the character to indicate alternative solutions.}
  \item{...}{additional arguments accepted by the function \code{\link{colorgram}} for \code{cp.plot}, or \code{\link{contour}} for \code{cp.plot.contours}.}
  \item{level}{decimal value of the confidence interval}
}

\value{
  The \code{get.ci} function returns a list with 6 items:
  \item{values.ploidy}{matrix of ploidy values with respective posterior probability.}
  \item{confint.ploidy}{boundaries of the confidence interval of the estimated ploidy.}
  \item{max.ploidy}{point estimate of the ploidy value that has the maximum posterior probability.}
  \item{values.cellularity}{matrix of cellularity values with respective posterior probability.}
  \item{confint.cellularity}{boundaries of the confidence interval of the estimated cellularity.}
  \item{max.cellularity}{point estimate of the cellularity value that has the maximum posterior probability.}
}

\examples{

data(CP.example)
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)

# Plot more contours
cp.plot(CP.example)
cp.plot.contours(CP.example, likThresh = c(0.95, 0.9999), add = TRUE)

# Return the 95% confidence interval
CP.example.ci <- get.ci(CP.example)
str(CP.example.ci)

}
