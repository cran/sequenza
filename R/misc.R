dt2 <- function(x, df, ncp, log = FALSE, mean, sd) {
    x2 <- (x - mean) / sd
    dt(x2, df = df, ncp = ncp, log = log)
}

split_chr_coord <- function (x) {
    # Ensure that there is a start and a end coordinate
    split_chr <- strsplit(x, split = ":")[[1]]
    chromosome <- split_chr[1]
    start_end <- split_chr[2]
    split_coors <- strsplit(start_end, split = "-")[[1]]
    start <- split_coors[1]
    end <- split_coors[2]
    if (is.na(start)) {
        start <- 1
    }
    if (is.na(end)){
        end <- 2147483647
    }
    paste0(chromosome, ":", start, "-", end)
}

weighted.median <- function(x, w, na.rm=TRUE, ties=NULL) {
    if (missing(w)) {
        w <- rep(1, length(x))
    }
    if (na.rm == TRUE) {
        keep <- !(is.na(x) | is.na(w));
        x <- x[keep]
        w <- w[keep]
    } else if (any(is.na(x))) {
        return(NA)
    }

    if (any(w < 0)) {
        stop("The weight vactor can only contains positive numbers")
    }
    n <- length(w)
    keep <- (w > 0)
    nkeep <- sum(keep)
    if (nkeep < n) {
        x <- x[keep]
        w <- w[keep]
        n <- nkeep
    }
    wInfs <- is.infinite(w)
    if (any(wInfs)) {
        x <- x[wInfs]
        n <- length(x)
        w <- rep(1, n)
    }

    if (n == 0) {
        return(NA)
    }
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    wcum <- cumsum(w)
    wsum <- wcum[n]
    wmid <- wsum / 2
    lows <- (wcum <= wmid)
    k  <- sum(lows)

    if (k == 0) {
        return(x[1])
    }
    if (k == n) {
        return(x[n])
    }

    wlow  <- wcum[k]
    whigh <- wsum - wlow
    if (whigh > wmid) {
        return(x[k + 1])
    }
    (wlow * x[k] + whigh * x[k + 1]) / wsum
}
