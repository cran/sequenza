\name{example.seqz}
\docType{data}
\alias{example.seqz}

\title{Example \dQuote{seqz} data}

\description{
    The \dQuote{seqz} file is produced by \command{sequenza-utils} and
    typically has the file extension \file{.seqz}. The data here is
    representative of a seqz file derived from an exome-sequenced tumor sample,
    such as could be obtained from TCGA.
}

\usage{
    data(example.seqz)
}

\format{
    A data frame with 53937 rows and 14 columns:
    \tabular{rll}{
        [,1] \tab chromosome    \tab Chromosome name \cr
        [,2] \tab position      \tab Base position \cr
        [,3] \tab base.ref      \tab Base in the reference genome\cr
        [,4] \tab depth.normal  \tab Read depth in the normal sample\cr
        [,5] \tab depth.tumor   \tab Read depth in the tumor sample\cr
        [,6] \tab depth.ratio   \tab Ratio of \code{depth.tumor} and \code{depth.normal}\cr
        [,7] \tab Af            \tab A-allele frequency in the tumor sample\cr
        [,8] \tab Bf            \tab B-allele frequency in the tumor sample,
                                     in heterozygous positions only\cr
        [,9] \tab zygosity.normal\tab Zygosity of the normal sample:
                                     "hom" for homozygous or "het" for heterozygous\cr
        [,10] \tab GC.percent   \tab \% GC content\cr
        [,11] \tab good.reads   \tab Number of reads from the tumor sample
                                     which pass the quality threshold \cr
        [,12] \tab AB.normal    \tab Base(s) found in the normal sample,
                                     sorted by allele frequency if more than one\cr
        [,13] \tab AB.tumor     \tab Base(s) found in the tumor sample \emph{but not} in
                                     the normal specimen, with their observed
                                     frequencies, separated by colons \cr
        [,14] \tab tumor.strand \tab Identical to \code{AB.tumor} but indicating, for
                                     each variant base, the fraction of reads
                                     oriented in the forward direction \cr
    }
}

\source{
This is derived from a TCGA specimen, but has been scrambled to anonymize the source.
The reference genome is hg19.  The GC content was calculated in 50-base windows.
}

\details{
    \code{example.seqz} can be loaded in the standard R way via
    \code{data(example.seqz)}, or it can be read from a text file using
    \code{\link{read.seqz}}. The former is useful for examples and testing,
    whereas the latter is representative of the standard workflow.
}

\keyword{datasets}
