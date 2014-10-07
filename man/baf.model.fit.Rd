\name{baf.model.fit}
\alias{baf.model.fit}
\alias{mufreq.model.fit}
\title{Model fitting using Bayesian inference}

\description{
  Computes the log-posterior probability distribution for the specified range of cellularity and ploidy parameters
}

\usage{
  mufreq.model.fit(cellularity = seq(0.3, 1, by = 0.01),
                   ploidy = seq(1, 7, by = 0.1),
                   mc.cores = getOption("mc.cores", 2L), ...)
  baf.model.fit(cellularity = seq(0.3, 1, by = 0.01),
                ploidy = seq(1, 7, by = 0.1),
                mc.cores = getOption("mc.cores", 2L), ...)
}

\arguments{
  \item{cellularity}{vector of cellularity parameter values to be tested.}
  \item{ploidy}{vector of ploidy parameter values to be tested.}
  \item{mc.cores}{number of cores to use, defined as in \code{\link{mclapply}}.}
  \item{...}{any argument accepted by \code{\link{mufreq.bayes}} or \code{\link{baf.bayes}}.}
}

\value{
  A list of three items:
  \item{ploidy}{tested values of the ploidy parameter}
  \item{cellularity}{tested values of the cellularity parameter}
  \item{lpp}{log-posterior probability of each pair of cellularity/ploidy parameters.}
}

\details{
  \code{baf.model.fit} uses the function \code{\link{baf.bayes}} to infer the log-posterior probability of the model fit using the possible combinations of cellularity and ploidy values provided in the arguments.
  Similarly \code{mufreq.model.fit} fits the mutation/depth ratio model using the function \code{\link{mufreq.bayes}}.
  \code{baf.model.fit} is generally used to infer cellularity and ploidy on segmented chromosomes, while \code{mufreq.model.fit} can also be used to estimate these parameters from mutations, and in addition can be useful to give a rough estimate of sub-clonal fractions. Be aware  that the mutation frequency is more sensitive to sampling bias, and is generally more noisy compared to the segmented B-allele frequency.
}

\seealso{
  \code{\link{cp.plot}} for visualization of the resulting object, and \code{\link{get.ci}} to extract confidence intervals.
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
seqz.data  <- read.seqz(data.file, chr.name = 1)

# Correct the coverage of the loaded chromosome:
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /
                           gc.vect[as.character(seqz.data$GC.percent)]
# Select the heterozygous positions
seqz.hom  <- seqz.data$zygosity.normal == 'hom'
seqz.het  <- seqz.data[!seqz.hom, ]
# Detect breakpoints
breaks <- find.breaks(seqz.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
seg.s1 <- segment.breaks(seqz.data, breaks = breaks)

# filter out small ambiguous segments, and conveniently weight the segments by size:
seg.filtered <- seg.s1[(seg.s1$end.pos - seg.s1$start.pos) > 10e6, ]
weights.seg  <- 150 + round((seg.filtered$end.pos -
                             seg.filtered$start.pos) / 1e6, 0)
# get the genome wide mean of the normalized depth ratio:
avg.depth.ratio <- mean(gc.stats$adj[,2])
# run the BAF model fit

CP <- baf.model.fit(Bf = seg.filtered$Bf, depth.ratio = seg.filtered$depth.ratio,
                    weight.ratio = weights.seg,
                    weight.Bf = weights.seg,
                    avg.depth.ratio = avg.depth.ratio,
                    cellularity = seq(0.1,1,0.01),
                    ploidy = seq(0.5,3,0.05))

confint <- get.ci(CP)
ploidy   <- confint$max.ploidy
cellularity <- confint$max.cellularity

   }
}
