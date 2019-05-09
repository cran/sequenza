context("Test read seqz data")
seqz_file <- system.file("extdata", "example.seqz.txt.gz", package = "sequenza")
col_names <- c("chromosome", "position", "base.ref", "depth.normal",
              "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
              "GC.percent", "good.reads", "AB.normal", "AB.tumor",
              "tumor.strand")

test_that("Testing read seqz data", {
    full_file <- read.seqz(seqz_file)
    expect_equal(colnames(full_file), col_names)
})

test_that("Testing tbi presence", {
    tbi <- file.exists(paste(seqz_file, "tbi", sep = "."))
    expect_equal(tbi, TRUE)
})

test_that("Reading single chromosome", {
    chr_17 <- read.seqz(seqz_file, chr_name = "17")
    expect_equal(dim(chr_17), c(2988, 14))
})

test_that("Reading selected region", {
    chr_17_p <- read.seqz(seqz_file, chr_name = "17:7500000-7600000")
    expect_equal(dim(chr_17_p), c(13, 14))
})
