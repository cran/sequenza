\name{abf.data.abfreq}
\alias{abf.data.abfreq}

\title{ABfreq file example.}

\description{
  Example of ABfreq file: This data is representative of an exome-sequenced tumor sample, e.g. as could be obtained from TCGA.
}

\value{
ABfreq is a tab separated text file with column headers. The file has currently 13 columns:
\item{chromosome}{with the chromosome names}
\item{n.base}{with the base positions}
\item{base.ref}{with the base in the reference genome used (usually hg19). Note the \code{base.ref} is NOT the base of the germline.}
\item{depth.normal}{read depth observed in the normal sample}
\item{depth.sample}{read depth observed in the tumor sample}
\item{depth.ratio}{ratio of \code{depth.sample} and \code{depth.normal}}
\item{Af}{A-allele frequency observed in the tumor sample}
\item{Bf}{B-allele frequency observed in the tumor sample}
\item{ref.zygosity}{zygosity of the reference sample. "hom" corresponds to AA or BB, whereas "het" corresponds to AB or BA}
\item{GC.percent}{GC-content (percent), calculated from the reference genome in fixed nucleotide windows }
\item{good.s.reads}{number of reads that passed the quality threshold (threshold specified in the pre-processing software)}
\item{AB.germline}{base found in the germline sample}
\item{AB.sample}{base found in the tumor sample}
}

\usage{
  data(abf.data.abfreq)
}
