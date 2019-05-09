context("Test sequenza models")


test_that("Testing theoretical depth ratio", {
    r1 <- theoretical.depth.ratio(CNt = 3, cellularity = 0.5, ploidy = 2)
    r2 <- theoretical.depth.ratio(CNt = 3, cellularity = 0.8, ploidy = 2)
    r3 <- theoretical.depth.ratio(CNt = 4, cellularity = 0.8, ploidy = 4)
    expect_equal(r1, 1.25)
    expect_equal(r2, 1.4)
    expect_equal(r3, 1)
})

test_that("Testing theoretical BAF", {
    r1 <- theoretical.baf(CNt = 4, B = 2, cellularity = 0.5)
    r2 <- theoretical.baf(CNt = 3, B = 1, cellularity = 0.8)
    r3 <- theoretical.baf(CNt = 4, B = 0, cellularity = 0.8)
    expect_equal(r1, 0.5)
    expect_equal(round(r2, 3), 0.357)
    expect_equal(round(r3, 3), 0.056)
})

test_that("Testing theoretical mufreq", {
    r1 <- theoretical.mufreq(CNt = 4, Mt = 2, cellularity = 0.5)
    r2 <- theoretical.mufreq(CNt = 3, Mt = 1, cellularity = 0.8)
    r3 <- theoretical.mufreq(CNt = 4, Mt = 0, cellularity = 0.8)
    expect_equal(round(r1, 3), 0.333)
    expect_equal(round(r2, 3), 0.286)
    expect_equal(r3, 0)
})

test_that("Testing baf types matrix", {
    x <- 3
    CNts <- do.call(c, sapply(1:x, FUN = function(x) rep(x, x)))
    Bs <- do.call(c, sapply(1:x, FUN = function(x) {
            y <- x:1
            z <-  rep(x, x)
            z - y
        }
    ))
    mat_b <- data.frame(CNn = rep(2, length(Bs) + 1),
                   CNt = c(0, CNts), B = c(0, Bs))
    mat_b <- mat_b[mat_b[, "B"] <= mat_b[, "CNt"] / 2, ]
    r1 <- baf.types.matrix(CNt.min = 0, CNt.max = 3, CNn = 2)
    expect_equal(r1, mat_b)
})

test_that("Testing mufreq types matrix", {
    x <- 3
    CNts <- do.call(c, sapply(1:x, FUN = function(x) rep(x, x + 1)))
    Mts <- do.call(c, sapply(1:x, FUN = function(x) 0:x))
    mat_m <- data.frame(CNn = rep(2, length(Mts) + 1),
                   CNt = c(0, CNts), Mt = c(0, Mts))
    r1 <- mufreq.types.matrix(CNt.min = 0, CNt.max = 3, CNn = 2)
    expect_equal(r1, mat_m)
})
