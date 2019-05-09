\name{read.seqz}
\alias{read.seqz}
\alias{read.acgt}
\title{Read a seqz or acgt format file}

\description{
    Efficiently reads a seqz file into R.
}

\usage{
    read.seqz(file, n_lines = NULL, col_types = "ciciidddcddccc", chr_name = NULL,
        buffer = 33554432, parallel = 1,
        col_names = c("chromosome", "position", "base.ref", "depth.normal",
            "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
            "GC.percent", "good.reads", "AB.normal", "AB.tumor",
            "tumor.strand"),...)
}

\arguments{
    \item{file}{file name}
    \item{col_types}{a string describing the classes of each columns of the
        input file (see \code{\link{read_tsv}}). The default value corresponds
        to the columns of a seqz file.}
    \item{chr_name}{if specified, only the selected chromosome will be
        extracted instead of the entire file. For \code{tabix}-indexed files
        this argument can also be used to extract coordinated-selected genomic
        regions. E.g. \code{chr_name="5:1-1000000"} will select the first
        megabase of chromosome 5.}
    \item{n_lines}{vector of length 2 specifying the first and last line to
        read from the file. If specified, only the selected portion of the
        file will be used.}
    \item{buffer}{maximal size of each chunk in bytes(see
        \code{\link{chunk.apply}}).}
    \item{parallel}{integer, number of threads used to process a seqz file
        (see \code{\link{chunk.apply}}).}
    \item{col_names}{names of the columns of the seqz file. The default
        corresponds to the column names of a seqz file.}
    \item{...}{any arguments accepted by \code{read_tsv}.}
}

\details{
    \code{read.seqz} is a function that allows to efficiently access a
    \code{seqz} file by chromosome or by line numbers. The function can also
    access coordinate specific regions with \code{tabix}-indexed \code{seqz}
    files.
    The specific content of a \code{seqz} file is explained in the \code{value}
    section.
}

\format{
    A seqz file is a tab-separated text file with 14 columns and a header row.
    The first 3 columns are derived from the original \code{pileup}
    file and contain:
    \describe{
        \item{chromosome}{the chromosome name}
        \item{position}{the base position}
        \item{base.ref}{the base in the reference genome.
            Note that this is NOT necessarily the same base as in
            the normal specimen.}
    The remaining 10 columns contain the following information:
        \item{depth.normal}{read depth observed in the normal sample}
        \item{depth.tumor}{read depth observed in the tumor sample}
        \item{depth.ratio}{ratio of \code{depth.tumor} and \code{depth.normal}}
        \item{Af}{A-allele frequency observed in the tumor sample}
        \item{Bf}{B-allele frequency observed in the tumor sample in
            heterozygous positions}
        \item{zygosity.normal}{zygosity of the reference sample. "hom"
            corresponds to AA or BB, whereas "het" corresponds to AB or BA}
        \item{GC.percent}{GC-content (percent), calculated from the reference
            genome in fixed nucleotide windows }
        \item{good.reads}{number of reads that passed the quality threshold
            (threshold specified in the pre-processing software), in the
            tumor specimen}
        \item{AB.normal}{base(s) found in the germline sample; for heterozygous
            positions AB are sorted using the values of Af and Bf respectively}
        \item{AB.tumor}{base(s) found in the tumor sample not present in the
            normal specimen. The field include all the variants found in the
            tumor alignment, separated by a colon. Each variant contains the
            base and the observed frequency}
        \item{tumor.strand}{frequency of the variant nucleotides detected on
            the forward orientation. The field have a consistent structure with
            \code{AB.tumor}, indicating the fraction, relative to the total
            number of reads presenting the specific variant, orientated
            in the forward direction}
    }
}

\seealso{
  \code{read_delim}.
}

\examples{
   \dontrun{

    data_file <-  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")

    ## read chromosome 1 from an seqz file.
    seqz_data <- read.seqz(data_file, chr_name = 1)

    ## Fast access to chromosome X using the file metrics
    gc.stats <- gc.sample.stats(data_file)
    chrX <- gc.stats$file.metrics[gc.stats$file.metrics$chr == "X", ]
    seqz.data <- read.seqz(data_file, n_lines = c(chrX$start, chrX$end))

    ## Compare the running time of the two different methods.
    system.time(seqz.data <- read.seqz(data_file, n_lines = c(chrX$start, chrX$end)))
    system.time(seqz.data <- read.seqz(data_file, chr_name = "X"))

   }
}
