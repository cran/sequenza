get.ci <- function(cp.table, level = 0.95) {
    znormsort <- sort(cp.table$lpp, decreasing = TRUE)
    znormcumLik <- cumsum(znormsort)
    n <- sapply(level, function(x) sum(znormcumLik < x) + 1)
    LikThresh <- znormsort[n]
    values.x <- data.frame(x = cp.table$ploidy,
        y = apply(cp.table$lpp, 1, max))
    values.y <- data.frame(x = apply(cp.table$lpp, 2, max),
    y = cp.table$cellularity)
    up.x  <- max(values.x$x[values.x$y >= LikThresh])
    low.x <- min(values.x$x[values.x$y >= LikThresh])
    max.x <- values.x$x[which.max(values.x$y)]
    up.y  <- max(values.y$y[values.y$x >= LikThresh])
    low.y <- min(values.y$y[values.y$x >= LikThresh])
    max.y <- values.y$y[which.max(values.y$x)]
    values.x$y <- values.x$y / sum(values.x$y)
    values.y$x <- values.y$x / sum(values.y$x)
    results <- list()
    results$values.ploidy <- values.x
    results$confint.ploidy <- c(low.x, up.x)
    results$max.ploidy <- max.x
    results$values.cellularity <- values.y
    results$confint.cellularity <- c(low.y, up.y)
    results$max.cellularity <- max.y
    results
}

alternative.cp.solutions <- function(cp.table) {
    ci <- get.ci(cp.table)
    p.alt <- which(diff(sign(diff(ci$values.ploidy$y))) == -2) + 1
    get.alt <- function(idx.p, cp.table) {
        idx.c <- which.max(cp.table$lpp[idx.p, ])
        c(cellularity = cp.table$cellularity[idx.c],
            ploidy = cp.table$ploidy[idx.p],
            SLPP = cp.table$lpp[idx.p, idx.c])
        }
    res <- lapply(p.alt, FUN = function (x) get.alt(x, cp.table))
    res <- as.data.frame(do.call(rbind, res))
    if (nrow(res) > 0 ){
        res[order(res$SLPP, decreasing = TRUE), ]
    } else {
        data.frame(cellularity = ci$max.cellularity,
            ploidy = ci$max.ploidy,
            SLPP =  cp.table$lpp[which(cp.table$ploidy == ci$max.ploidy),
                which(cp.table$cellularity == ci$max.cellularity)])
   }
}
