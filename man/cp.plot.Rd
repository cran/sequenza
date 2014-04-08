\name{cp.plot}
\alias{cp.plot}
\alias{cp.plot.contours}
\title{Plot log-likelihood for the tested values of cellularity and ploidy}

\description{
  This function uses the \code{\link{colorgram}} function from the package \pkg{squash} to plot log-likelihood for the tested combinations of cellularity and ploidy
}

\usage{
   cp.plot(cp.table, xlab = "Ploidy", ylab = "Cellularity",
           zlab = "Scaled rank likelihood", 
           colFn = colorRampPalette(c('white', 'lightblue')), ...)
   cp.plot.contours(cp.table, likThresh = c(0.95),
                    col = palette(), legend.pos = "bottomright", pch = 18, ...)
}

\arguments{
  \item{cp.table}{list, as output from \code{\link{baf.model.fit}} or \code{\link{mufreq.model.fit}}.}
  \item{xlab}{xlab parameter as in the function \code{\link{colorgram}}.}
  \item{ylab}{ylab parameter as in the function \code{\link{colorgram}}.}
  \item{zlab}{zlab parameter as in the function \code{\link{colorgram}}.}
  \item{colFn}{colFn parameter as in the function \code{\link{colorgram}}.}  
  \item{likThresh}{vector of selected likelihood tresholds.}
  \item{col}{vector of colours.}
  \item{legend.pos}{position for placing the legend.}
  \item{pch}{character used to indicate the point estimate.}
  \item{...}{additional arguments accepted by the function \code{\link{colorgram}} for \code{cp.plot}, or \code{\link{contour}} for \code{cp.plot.contours}.}
}

\examples{

data(CP.example)
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)

data(CPsnp.example)
cp.plot(CPsnp.example)
cp.plot.contours(CPsnp.example, likThresh = c(0.95, 0.9999), add = TRUE)

}
