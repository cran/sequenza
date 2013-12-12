\name{sequenza}
\alias{sequenza.extract}
\alias{sequenza.fit}
\title{Use sequenza to estimate tumor purity and ploidy.}

\description{
  The main interface of the package, to run several of the functions of \code{sequenza} in a standardized pipeline.
}

\usage{
  sequenza.extract(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80,
                   kmin = 10, mufreq.treshold = 0.10, min.reads = 40,
                   max.mut.types = 1, min.type.freq = 0.9)

  sequenza.fit(sequenza.extract, female = TRUE, segment.filter = 1e7,
               XY = c(X = "X", Y = "Y"), cellularity = seq(0.1,1,0.01),
               ploidy = seq(1, 7, 0.1), priors.table = data.frame(CN = 2, value = 2),
               mc.cores = getOption("mc.cores", 2L))
}

\arguments{
  \item{file}{an ABfreq file.}
  \item{gz}{logical. If TRUE (the default) the function expects a gzipped file.}
  \item{window}{size of windows used when plotting mean and quartile ranges of depth ratios and B-allele frequencies. Smaller windows will take more time to compute.}
  \item{overlap}{integer specifying the number of overlapping windows.}
  \item{gamma, kmin}{arguments passed to \code{aspcf} from the \code{copynumber} package.}
  \item{mufreq.treshold}{mutation frequency threshold.}
  \item{min.reads}{minimal number of reads above the quality threshold to accept the mutation call.}
  \item{max.mut.types}{maximal number of different base substitutions per position. Integer from 1 to 3 (since there are only 4 bases). Default is 3, to accept "noisy" mutation calls.}
  \item{min.type.freq}{minimal frequency of aberrant types.}
  \item{sequenza.extract}{a list of objects as output from the \code{sequenza.extract} function.}
  \item{female}{logical, indicating whether the sample is male or female, to properly handle the X and Y chromosomes. Implementation only works for the human normal karyotype.}
  \item{segment.filter}{threshold segment length (in base pairs) to filter out short segments, that can cause noise when fitting the cellularity and ploidy parameters. The threshold will not affect the allele-specific segmentation.}
  \item{XY}{character vector of length 2 specifying the labels used for the X and Y chromosomes. Defaults to c(X = "X", Y = "Y").}
  \item{cellularity}{vector of values to test as cellularity parameter.}
  \item{ploidy}{vector of values to test as ploidy parameter.}
  \item{priors.table}{data frame with the columns \code{CN} and \code{value}, containing the copy numbers and the corresponding weights. To every copy number is assigned the value 1 as default, so every values different then 1 will change the corresponding weight.}
  \item{mc.cores}{number of cores to use, defined as in the \code{parallel} package.}
}

\details{
   The function \code{sequenza.extract} utilizes a range of functions from the sequenza package to read the raw data, normalize the depth.ratio for GC-content bias, perform allele-specific segmentation, filter for noisy mutations and binning of the raw data for plotting. The computed objects are returned as a single list object. This object can be given to  \code{sequenza.fit}, which uses \code{baf.model.fit} to calculate the log-likelihood for all pairs of the ploidy and cellularity parameters.
}

\seealso{
  \code{genome.view}, \code{baf.bayes}, \code{cp.plot}, \code{get.ci}.
}

\examples{
  \dontrun{

data.file <-  system.file("data", "abf.data.abfreq.txt.gz",
              package = "sequenza")
test <- sequenza.extract(data.file)
CP   <- sequenza.fit(test, mc.cores = 4)
   }
}
