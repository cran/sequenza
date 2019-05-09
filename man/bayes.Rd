\name{baf.bayes}
\alias{baf.bayes}
\alias{mufreq.bayes}
\title{Model allele-specific copy numbers with specified cellularity and ploidy parameters}

\description{
  Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele frequency and depth ratio.
}

\usage{
  baf.bayes(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            sd.Bf = 0.1, sd.ratio = 0.5, weight.Bf = 1, weight.ratio = 1,
            CNt.min = 0, CNt.max = 7, CNn = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1),
            ratio.priority = FALSE)
  mufreq.bayes(mufreq, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7, CNn = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1))
}

\arguments{
  \item{Bf}{vector of B-allele frequencies (values can range from 0 to 0.5).}
  \item{mufreq}{vector of mutation frequencies (values can range from 0 to 1).}
  \item{depth.ratio}{vector of depth ratios.}
  \item{sd.ratio}{standard deviation observed in the depth ratio measures in a segment}
  \item{sd.Bf}{standard deviation observed in the B-allele frequency measures in a segment}
  \item{weight.Bf}{vector of weights for B-allele frequency values.}
  \item{weight.mufreq}{vector of weights for the mutation frequency values.}
  \item{weight.ratio}{vector of weights for the depth ratio values.}
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{2 * ratio between total DNA content in a tumor cell and a normal cell.}
  \item{avg.depth.ratio}{average normalized depth ratio.}
  \item{CNt.min}{minimum copy number to consider in the model.}
  \item{CNt.max}{maximum copy number to consider in the model.}
  \item{CNn}{copy number of the normal genome.}
  \item{priors.table}{data frame with columns \code{CN} and \code{value}, containing the copy numbers and the corresponding weights. To every copy number is assigned the value 1 as default, so any values different from 1 will change the corresponding weight.}
  \item{ratio.priority}{logical, if TRUE only the depth ratio will be used to determine the copy number state, while the Bf value will be used to determine the number of B-alleles.}
}

\value{
  \item{CNt}{copy number of the tumor cell at the tested point.}
  \item{A}{number of A-alleles at the tested point.}
  \item{B}{number of B-alleles at the tested point.}
  \item{CNn}{copy number of the normal cell at the tested point (equal to CNn given as argument).}
  \item{Mt}{number of mutated alleles at the tested point.}
  \item{LPP}{log-posterior probability of model fitting at the given point/segment.}
}
\details{
  \code{baf.bayes} and \code{mufreq.bayes} use a naive Bayesian approach to calculate the posterior probability of fitness of the data point with the model point resulting from the given values of cellularity and DNA-content.
}

\seealso{
   \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
  \dontrun{
data.file <-  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")
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

#detect copy number alteration on the segments:

cn.alleles <- baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                        cellularity = cellularity, ploidy = ploidy,
                        avg.depth.ratio = 1)

head(cbind(seg.s1, cn.alleles))

# create mutation table:
mut.tab   <- mutation.table(seqz.data, mufreq.treshold = 0.15,
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9, segments = seg.s1)

mut.tab.clean <- na.exclude(mut.tab)

# Detect mutated alleles:
mut.alleles <- mufreq.bayes(mufreq = mut.tab.clean$F,
                            depth.ratio = mut.tab.clean$adjusted.ratio,
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio)
head(cbind(mut.tab.clean[,c("chromosome","position","F",
                            "adjusted.ratio", "mutation")],
           mut.alleles))


   }
}