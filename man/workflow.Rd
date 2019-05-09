\name{sequenza}
\alias{sequenza.extract}
\alias{sequenza.fit}
\alias{sequenza.results}
\title{Sequenza convenience functions for standard analysis}

\description{
  These three functions are intended to be the main user interface of the package, to run several of the functions of \code{sequenza} in a standardized pipeline.
}

\usage{
  sequenza.extract(file, window = 1e6, overlap = 1,
    gamma = 80, kmin = 10, gamma.pcf = 140, kmin.pcf = 40,
    mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
    breaks = NULL, breaks.method = "het", assembly = "hg19",
    weighted.mean = TRUE, normalization.method = "mean",
    ignore.normal = FALSE, parallel = 1, gc.stats = NULL,
    segments.samples = FALSE)

  sequenza.fit(sequenza.extract, female = TRUE, N.ratio.filter = 10,
               N.BAF.filter = 1, segment.filter = 3e6,
               mufreq.treshold = 0.10, XY = c(X = "X", Y = "Y"),
               cellularity = seq(0.1,1,0.01), ploidy = seq(1, 7, 0.1),
               ratio.priority = FALSE, method = "baf",
               priors.table = data.frame(CN = 2, value = 2),
               chromosome.list = 1:24, mc.cores = getOption("mc.cores", 2L))

  sequenza.results(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                   cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                   ratio.priority = FALSE, XY = c(X = "X", Y = "Y"),
                   chromosome.list = 1:24)
}

\arguments{
  \item{file}{the name of the seqz file to read.}
  \item{window}{size of windows used when plotting mean and quartile ranges of depth ratios and B-allele frequencies. Smaller windows will take more time to compute.}
  \item{overlap}{integer specifying the number of overlapping windows.}
  \item{gamma, kmin}{arguments passed to \code{\link{aspcf}} from the \pkg{copynumber} package.}
  \item{gamma.pcf, kmin.pcf}{arguments passed to \code{\link{pcf}} from the \pkg{copynumber} package. The arguments are effective only when \code{breaks.method} is set to "full".}
  \item{mufreq.treshold}{mutation frequency threshold.}
  \item{min.reads}{minimum number of reads above the quality threshold to accept the mutation call.}
  \item{min.reads.normal}{minimum number of reads used to determine the genotype in the normal sample.}
  \item{min.reads.baf}{threshold on the depth of the positions included to calculate the average BAF for segment.}
  \item{max.mut.types}{maximum number of different base substitutions per position. Integer from 1 to 3 (since there are only 4 bases). Default is 3, to accept "noisy" mutation calls.}
  \item{min.type.freq}{minimum frequency of aberrant types.}
    \item{min.fw.freq}{minimum frequency of variant reads detected in the forward strand. Setting it to 0, all the variant calls with strand frequency in the interval outside 0 and 1, margin not comprised, would be discarded.}
  \item{verbose}{logical, indicating whether to print information about the chromosome being processed.}
  \item{chromosome.list}{vector containing the index or the names of the chromosome to include in the model fitting.}
  \item{breaks}{Optional data.frame in the format chrom, start.pos, end.pos, defining a pre-existing segmentation. When the argument is set the built-in segmentation will be skipped in favor of the suggested breaks.}
  \item{breaks.method}{Argument indicating the resolution of the segmentation. Possible values are \code{fast}, \code{het} and \code{full}, where \code{fast} allows the lower resolution and \code{full} the higher. Custom values of \code{gamma} and \code{kmin} need to be adjusted to have optimal results.}
  \item{assembly}{assembly version of the genome, see \code{\link{aspcf}} or \code{\link{pcf}}.}
  \item{weighted.mean}{boolean to select if the segments should be calculated using the read depth as weights to calculate depth ratio and B-allele frequency means.}
  \item{normalization.method}{string defining the operation to perform during the GC-normalization process. Possible values are \code{mean} (default) and \code{median}. A \code{median} normalization is preferable with noisy data.}
  \item{ignore.normal}{boolean, when set to TRUE the process will ignore the normal coverage and perform the analysis by using the normalized tumor coverage.}
  \item{parallel}{integer, number of threads used to process a seqz file
        (see \code{\link{chunk.apply}}).}
  \item{gc.stats}{object returned from the function \code{\link{gc.sample.stats}}. If \code{NULL} the object will be computed from the input file.}
  \item{segments.samples}{EXPERIMENTAL. Segment both tumor and normal samples separately, and add it to the QC plots.}
  \item{sequenza.extract}{a list of objects as output from the \code{sequenza.extract} function.}
  \item{method}{method to use to fit the data; possible values are \code{baf} to use \code{\link{baf.model.fit}} or \code{mufreq} to use the \code{\link{mufreq.model.fit}} function to fit the data.}
  \item{cp.table}{a list of objects as output from the \code{sequenza.fit} function.}
  \item{female}{logical, indicating whether the sample is male or female, to properly handle the X and Y chromosomes. Implementation only works for the human normal karyotype.}
  \item{CNt.max}{maximum copy number to consider in the model.}
  \item{N.ratio.filter}{threshold of minimum number of observation of depth ratio in a segment.}
  \item{N.BAF.filter}{threshold of minimum number of observation of B-allele frequency in a segment.}
  \item{segment.filter}{threshold segment length (in base pairs) to filter out short segments, that can cause noise when fitting the cellularity and ploidy parameters. The threshold will not affect the allele-specific segmentation.}
  \item{XY}{character vector of length 2 specifying the labels used for the X and Y chromosomes.}
  \item{cellularity}{vector of candidate cellularity parameters.}
  \item{ploidy}{vector candidate ploidy parameters.}
  \item{priors.table}{data frame with the columns \code{CN} and \code{value}, containing the copy numbers and the corresponding weights. To every copy number is assigned the value 1 as default, so every values different then 1 will change the corresponding weight.}
  \item{ratio.priority}{logical, if TRUE only the depth ratio will be used to determine the copy number state, while the Bf value will be used to determine the number of B-alleles.}
  \item{sample.id}{identifier of the sample, to be used as a prefix for saved objects.}
  \item{out.dir}{output directory where the files and objects will be saved.}
  \item{mc.cores}{legacy argument to set the number of cores, but it refers to the \code{cl} of \code{\link{pblapply}}. It uses \code{\link{mclapply}} when set to an integer.}

}

\details{
   The first function, \code{sequenza.extract}, utilizes a range of functions from the sequenza package to read the raw data, normalize the depth.ratio for GC-content bias, perform allele-specific segmentation, filter for noisy mutations and bin the raw data for plotting. The computed objects are returned as a single list object.

   The segmentation by default is performed using only the heterozygous position and the \code{\link{aspcf}} function from \pkg{copynumber} package. The \code{full} option in the \code{breaks.method} argument allow to combine results of the segmentation of all the data available, using the \code{\link{pcf}} function, and the default \code{\link{aspcf}} using only the heterozygous positions.

   The second function, \code{sequenza.fit}, accepts the output from \code{sequenza.extract} and calls \code{\link{baf.model.fit}} to calculate the log-posterior probability for all pairs of the candidate ploidy and cellularity parameters.

   The third function, \code{sequenza.results}, saves a number of objects in a specified directory (default is the working directory). The objects are:
\itemize{
 \item The list of segments with resulting copy numbers and major and minor alleles.
 \item The candidate mutation list with variant allele frequency, and copy number and number of mutated allele, in relation of the clonal population (for sub-clonal population it needs to be processed with further methods).
 \item A plot of all the chromosomes in one image, representing the major and minor alleles and the absolute copy number changes (genome_view).
 \item Multiple plots with one chromosome per image, representing copy-number, B-allele frequency and mutation in parallel (chromosome_view).
 \item Results of the model fitting (CP_contours and confints_CP)
 \item A summary of the copy number state of the sample (CN_bars).
}
}

\seealso{
  \code{\link{genome.view}}, \code{\link{baf.bayes}}, \code{\link{cp.plot}}, \code{\link{get.ci}}.
}

\examples{
  \dontrun{

data.file <-  system.file("extdata", "example.seqz.txt.gz",
              package = "sequenza")
test <- sequenza.extract(data.file)
test.CP   <- sequenza.fit(test)
sequenza.results(test, test.CP, out.dir = "example",
                 sample.id = "example")

   }
}