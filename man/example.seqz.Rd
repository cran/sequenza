\name{example.seqz}
\docType{data}
\alias{example.seqz}

\title{Example \dQuote{seqz} data}

\description{
 The \dQuote{seqz} file is produced by \command{sequenza-utils.py} and typically has the file extension \file{.seqz}. The data here is representative of an exome-sequenced tumor sample, such as could be obtained from TCGA. 
}

\usage{
  data(example.seqz)
}

\details{
\code{example.seqz} can be loaded in the standard R way via \code{data(example.seqz)}, or it can be read from a text file using \code{\link{read.seqz}}. The former is useful for examples and testing, whereas the latter is representative of the standard workflow.
}

\format{
A data frame with a header row and 13 columns:
\describe{
\item{chromosome}{with the chromosome name}
\item{position}{with the base position}
\item{base.ref}{with the base in the reference genome used (usually hg19). Note the \code{base.ref} is NOT necessarily the base in the normal specimen.}
The remaining 10 columns contain the following information:
\item{depth.normal}{read depth observed in the normal sample}
\item{depth.tumor}{read depth observed in the tumor sample}
\item{depth.ratio}{ratio of \code{depth.tumor} and \code{depth.normal}}
\item{Af}{A-allele frequency observed in the tumor sample}
\item{Bf}{B-allele frequency observed in the tumor sample in heterozygous positions}
\item{zygosity.normal}{zygosity of the reference sample. "hom" corresponds to AA or BB, whereas "het" corresponds to AB or BA}
\item{GC.percent}{GC-content (percent), calculated from the reference genome in fixed nucleotide windows }
\item{good.reads}{number of reads that passed the quality threshold (threshold specified in the pre-processing software), in the tumor specimen}
\item{AB.normal}{base(s) found in the germline sample; for heterozygous positions AB are sorted using the values of Af and Bf respectively}
\item{AB.tumor}{base(s) found in the tumor sample not present in the normal specimen. The field include all the variants found in the tumor alignment, separated by a colon. Each variant contains the base and the observed frequency}
\item{tumor.strand}{frequency of the variant nucleotides detected on the forward orientation. The field have a consistent structure with \code{AB.tumor}, indicating the fraction, relative to the total number of reads presenting the specific variant, orientated in the forward direction}
}

}

\keyword{datasets}
