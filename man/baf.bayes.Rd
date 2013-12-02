\name{baf.bayes}
\alias{baf.bayes}
\alias{mufreq.bayes}
\title{Model allele-specific copy numbers with specified cellularity and DNA-content parameters}

\description{
  Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-likelihood of the fit, for given values of B-allele frequency and depth ratio.
}

\usage{
  baf.bayes(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            weight.Bf = 100, weight.ratio = 100, CNt.min = 0, CNt.max = 7,
            CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max, value = 1),
            ratio.priority = FALSE, skew.baf = 0.95)
  mufreq.bayes(mufreq, depth.ratio, cellularity, ploidy, avg.depth.ratio,
            weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7, CNn = 2,
            priors.table = data.frame(CN = CNt.min:CNt.max, value = 1))
}

\arguments{
  \item{Bf}{vector of B-allele frequencies (values can range from 0 to 0.5).}
  \item{mufreq}{vector of mutation frequencies (values can range from 0 to 1).}
  \item{depth.ratio}{vector of depth ratios.}
  \item{weight.Bf}{vector of weights for B-allele frequency values.}
  \item{weight.mufreq}{vector of weights for the mutation frequency values.}  
  \item{weight.ratio}{vector of weights for the depth ratio values.}
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{ratio between total DNA content in a tumor cell and a normal cell.}
  \item{avg.depth.ratio}{average normalized depth ratios.}
  \item{CNt.min}{minimum copy number to consider in the model.}
  \item{CNt.max}{maximum copy number to consider in the model.}
  \item{CNn}{copy number of the normal genome.}
  \item{priors.table}{data frame with the columns \code{CN} and \code{value}, containing the copy numbers and the corresponding weights. To every copy number is assigned the value 1 as default, so every values different then 1 will change the corresponding weight.}
  \item{ratio.priority}{logical, if TRUE only the depth ratio will be used to determine the copy number state, while the Bf value will be used to determine the number of B-alleles.}
  \item{skew.baf}{observed B-allele frequency often does not converge to the 0.5 frequency, for normal diploid positions. This argument indicate at which percentile of the observed B-allele frequency to skew the theoretic value below 0.5.}
}

\value{
  \item{CNt}{copy number of the tumor cell at the tested point.}
  \item{A}{number of A-alleles at the tested point.}
  \item{B}{number of B-alleles at the tested point.}
  \item{CNn}{copy number of the normal cell at the tested point (equal to CNn given as argument).}
  \item{Mt}{number of mutated alleles at the tested point.}
  \item{L}{log-likelihood of model fitting at the given point.}
}
\details{
  \code{baf.bayes} and \code{mufreq.bayes} use a naive Bayesian approach to calculate the likelihood of fitness of the data point with the model point resulting from the given values of cellularity and DNA-content.
}

\seealso{
   \code{baf.model.fit}, \code{mufreq.model.fit}.
}

\examples{
  \dontrun{
data.file <-  system.file("data", "abf.data.abfreq.txt.gz", package = "sequenza")
# read all the chromosomes:
abf.data  <- read.abfreq(data.file)
# Gather genome wide GC-stats from raw file:
gc.stats <- gc.sample.stats(data.file)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
# Read only one chromosome:
abf.data  <- read.abfreq(data.file, chr.name = 1)

# Correct the coverage of the loaded chromosome:
abf.data$adjusted.ratio <- abf.data$depth.ratio / 
                           gc.vect[as.character(abf.data$GC.percent)]
# Select the heterozygous positions
abf.hom  <- abf.data$ref.zygosity == 'hom'
abf.het  <- abf.data[!abf.hom, ]
# Detect breakpoints
breaks <- find.breaks(abf.het, gamma = 80, kmin = 10, baf.thres = c(0, 0.5))
# use heterozygous and homozygous position to measure segment values
seg.s1 <- segment.breaks(abf.data, breaks = breaks)

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
ploidy   <- confint$max.x
cellularity <- confint$max.y

#detect copy number alteration on the segments:

cn.alleles <- baf.bayes(Bf = seg.s1$Bf, depth.ratio = seg.s1$depth.ratio,
                        cellularity = cellularity, ploidy = ploidy,
                        avg.depth.ratio = 1)

head(cbind(seg.s1, cn.alleles))

# create mutation table:
mut.tab   <- mutation.table(abf.data, mufreq.treshold = 0.15, 
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9, segments = seg.s1)
                            
mut.tab.clean <- na.exclude(mut.tab)

# Detect mutated alleles:
mut.alleles <- mufreq.bayes(mufreq = mut.tab.clean$F,
                            depth.ratio = mut.tab.clean$adjusted.ratio,
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio)
head(cbind(mut.tab.clean[,c("chromosome","n.base","F","adjusted.ratio", "mutation")], mut.alleles))


   }
}