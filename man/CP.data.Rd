\name{CP.example}
\docType{data}
\alias{CP.example}

\title{Example of cellularity and ploidy results}

\description{
    Examples of results from the maximum a posteriori estimation from a set
    of cellularity and ploidy values, as returned by the functions
    \code{\link{baf.model.fit}} and \code{\link{mufreq.model.fit}}.
}

\usage{
    data(CP.example)
}

\format{
    A list containing three items:
        \describe{
        \item{ploidy}{numeric vector of tested ploidy values.}
        \item{cellularity}{numeric vector of tested cellularity values.}
        \item{lpp}{numeric matrix of log-posterior probability for each
        \var{(ploidy, cellularity)} pair.}
    }
}

\examples{
    data(CP.example)
    str(CP.example)

    ## Visualization of the object
    image(x = CP.example$ploidy,
        y = CP.example$cellularity,
        z = CP.example$lpp)

    ## A better plot
    cp.plot(CP.example)
    cp.plot.contours(CP.example, add = TRUE)
}

\keyword{datasets}
