## ----load----------------------------------------------------------------
library(sequenza)

## ----data----------------------------------------------------------------
data.file <-  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")

## ----cp_data, echo=FALSE, message=FALSE----------------------------------
library(sequenza)
data(CP.example)
CP <- CP.example

## ----extract, message=FALSE, warning=FALSE, results="hide"---------------
test <- sequenza.extract(data.file, verbose = FALSE)

## ----fit, eval=FALSE-----------------------------------------------------
#  CP <- sequenza.fit(test)

## ----results-------------------------------------------------------------
sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = "Test",
    out.dir="TEST")

## ----res_dir, echo=FALSE-------------------------------------------------

# create results list 
# 
res_list <- c("alternative_fit.pdf", "alternative_solutions.txt",
    "chromosome_depths.pdf", "chromosome_view.pdf",
    "CN_bars.pdf", "confints_CP.txt",
    "CP_contours.pdf", "gc_plots.pdf",
    "genome_view.pdf", "model_fit.pdf",
    "mutations.txt", "segments.txt",
    "sequenza_cp_table.RData", "sequenza_extract.RData",
    "sequenza_log.txt")
res_list <- paste("Test", res_list, sep = "_")

description_list <- c(
    "Alternative solution fir to the segments. One solution per slide",
    "List of all ploidy/cellularity alternative solution",
    "Visualization of sequencing coverage in the normal and in the tumor samples, before and after normalization",
    "Visualization per chromosome of depth.ratio, B-allele frequency and mutations, using the selected or estimated solution. One chromosome per slide",
    "Bar plot representing the percentage of genome in the detected copy number states",
    "Table of the confidence inerval of the best solution from the model",
    "Visualization of the likelihood density for each pair of cellularity/ploidy solution. The local maximum-likelihood points and confidence interval of the best estimate are also visualized",
    "Visualization of the GC correction in the normal and in the tumor sample",
    "Genome-whide visualization of the allele-specific and absolute copy number results, and raw profile of the depth ratio and allele frequency",
    "model_fit.pdf",
    "Table with mutation and estimated number of mutated alleles (Mt)",
    "Table listing the detected segments, with estimated copy number state at each sement",
    "RData object dump of the maxima a posteriori computation",
    "RData object dump of all the sample information",
    "Log with version and time information")
knitr::kable(data.frame(Files = res_list, Description = description_list))
#dir("TEST", pattern = "Test")


## ----read_segs, echo=FALSE-----------------------------------------------
seg.tab <- read.table("TEST/Test_segments.txt",
                      header = TRUE, sep ="\t")

alt_res <- read.table("TEST/Test_alternative_solutions.txt",
                      header = TRUE, sep ="\t")
seg.tab <- seg.tab[seg.tab$CNt <= 4, ]
is.num <- sapply(seg.tab, is.numeric)
seg.tab[is.num] <- lapply(seg.tab[is.num], round, 3)

## ----head_segs, echo=FALSE-----------------------------------------------
knitr::kable(head(seg.tab))


## ----g_view, echo=FALSE, fig.height=5, fig.width=10, fig.align='center'----
sequenza:::genome.view(seg.tab)

## ----g_view_tot, echo=FALSE, fig.height=5, fig.width=10, fig.align='center'----
sequenza:::genome.view(seg.tab, info.type = "CNt")

## ----g_view_raw, echo=FALSE, fig.height=5, fig.width=10, fig.align='center'----
sequenza:::plotRawGenome(test, cellularity = alt_res$cellularity[1],
    ploidy = alt_res$ploidy[1])

## ----CPplot, echo=TRUE, fig.height=5, fig.width=5, fig.align='center'----
cp.plot(CP)
cp.plot.contours(CP, add = TRUE,
   likThresh = c(0.999, 0.95),
   col = c("lightsalmon", "red"), pch = 20)

## ----c_view, echo=TRUE, fig.height=6, fig.width=8, fig.align='center'----
chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]],
                ratio.windows = test$ratio[[1]],  min.N.ratio = 1,
                segments = test$segments[[1]],
                main = test$chromosomes[1],
                cellularity = 0.89, ploidy = 1.9,
                avg.depth.ratio = 1)

