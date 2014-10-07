\name{chromosome.view}
\alias{chromosome.view}
\alias{genome.view}

\title{A graphical representation of multiple chromosomal features in parallel}

\description{
  A graphical representation of multiple chromosomal features in parallel.
}

\usage{
chromosome.view(baf.windows, ratio.windows, mut.tab = NULL, segments = NULL,
                min.N.baf = 1, min.N.ratio = 10000, main = "", vlines = FALSE,
                legend.inset = c(-20 * strwidth("a", units = "figure"), 0),
                BAF.style = "lines", CNn = 2, cellularity = NULL, ploidy = NULL,
                avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24",
                model.col = 1, x.chr.space = 10)
genome.view(seg.cn, info.type = "AB", ...)
}

\arguments{
  \item{baf.windows}{matrix containing the windowed B-allele frequency values for one chromosome.}
  \item{ratio.windows}{matrix containing the windowed depth ratio values for one chromosome.}
  \item{mut.tab}{mutation table of one chromosome. If specified, the mutations will be drawn in a top panel. \code{mut.tab} must be output from the \code{\link{mutation.table}} function.}
  \item{segments}{segmentation for one chromosome. If specified, the segmented B-allele frequency and depth ratio values will be shown as red lines.}
  \item{min.N.baf}{minimum number of observations required in a BAF window for plotting.}
  \item{min.N.ratio}{minimum number of observations required in a depth ratio window for plotting.}
  \item{CNn}{copy number of the germline genome. }
  \item{vlines}{logical, if TRUE the plot will include dotted vertical lines corresponding to segment breaks.}
  \item{cellularity}{fraction of tumor cells in the sample.}
  \item{ploidy}{value of the estimated \code{ploidy} parameter.}
  \item{avg.depth.ratio}{the average value of the normalized depth ratio.}
  \item{main}{main title of the plot.}
  \item{legend.inset}{the inset argument to pass to the \code{\link{legend}} function. Defines the distance between the mutation legend and the plot border.}
  \item{BAF.style}{the style to visualize the theoretical thresholds for the B allele frequency panel. Possible values are "lines", "blocks", "none". Default is "none".}
  \item{model.lwd}{width of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.lty}{line type of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{model.col}{color of the theoretical lines, if the segments matrix contains the columns A, B and CNt.}
  \item{x.chr.space}{step in megabase on the positions to visualize on the x-axis.}
  \item{seg.cn}{genome wide segments, with the columns A, B and CNt.}
  \item{info.type}{information to plot in \code{genome.view}. Available options are "CNt" to draw the copy numbers and "AB" (default) to draw the two alleles.}
  \item{...}{optional arguments passed to \code{\link{plot}}.}
}

\details{
  \code{chromosome.view} is a convenient plotting function using \code{\link{plotWindows}} and some handy \code{\link{par}} parameters, to display the raw sequencing results along with the model fitting.
}

\seealso{
  \code{\link{windowValues}}, \code{\link{find.breaks}}.
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

# Binning the values of depth ratio and B allele frequency
seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
             positions = seqz.data$position,
             chromosomes = seqz.data$chromosome,
             window = 1e6, overlap = 1,
             weight = seqz.data$depth.normal)

seqz.b.win <- windowValues(x = seqz.het$Bf,
             positions = seqz.het$position,
             chromosomes = seqz.het$chromosome,
             window = 1e6, overlap = 1,
             weight = round(x = seqz.het$good.reads, digits = 0))
# create mutation table:
mut.tab   <- mutation.table(seqz.data, mufreq.treshold = 0.15,
                            min.reads = 40, max.mut.types = 1,
                            min.type.freq = 0.9, segments = seg.s1)
# chromosome view without parametes:
chromosome.view(mut.tab = mut.tab[mut.tab$chromosome == "1",],
                baf.windows = seqz.b.win[[1]],
                ratio.windows = seqz.r.win[[1]], min.N.ratio = 1,
                segments = seg.s1[seg.s1$chromosome == "1",],
                main = "Chromosome 1")

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

seg.s1 <- cbind(seg.s1, cn.alleles)

# Chromosome view with estimated paramenters:
chromosome.view(mut.tab = mut.tab[mut.tab$chromosome == "1",],
                baf.windows = seqz.b.win[[1]],
                ratio.windows = seqz.r.win[[1]], min.N.ratio = 1,
                segments = seg.s1[seg.s1$chromosome == "1",],
                main = "Chromosome 1", cellularity = cellularity,
                ploidy = ploidy, avg.depth.ratio = 1,
                BAF.style = "lines")
   }
}
