windowValues <- function(x, positions, chromosomes, window = 1e6,
                         overlap = 0, weight = rep.int(x = 1,
                         times = length(x)), start.coord = 1) {
    weight   <- sqrt(weight)
    overlap  <- as.integer(overlap)
    window.offset <- as.integer(window -
        round(window * (overlap / (overlap + 1))))
    chr.ordered <- unique(chromosomes)
    data.splitByChr    <- split(data.frame(pos = positions, x = x,
                                           weight = weight),
                                f = factor(chromosomes, levels = chr.ordered))
    lapply(data.splitByChr, function(x) {
        range.pos <- range(x$pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2]) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        pos.cut <- cut(x$pos, breaks = beam.coords)
        x.split <- split(x$x, f = pos.cut)
        weight.split <- split(x$weight, f = pos.cut)
        window.starts <- beam.coords[1:nWindows]
        window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        x.window <- lapply(idx.list, function(idx) {
            unlist(x.split[idx], use.names = FALSE)
        })
        weight.window <- lapply(idx.list, function(idx) {
            unlist(weight.split[idx], use.names = FALSE)
        })
        window.means <- mapply(weighted.mean, x = x.window, w = weight.window)
        window.quantiles <- sapply(x.window, quantile,
            probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
        window.counts <- sapply(x.window, length)
        data.frame(start = window.starts, end = window.ends,
            mean = window.means, q0 = window.quantiles[1, ],
            q1 = window.quantiles[2, ], N = window.counts)
    })
}

windowBf <- function(Af, Bf, good.reads, positions, chromosomes,
    window = 1e6, overlap = 0, start.coord = 1, conf = 0.95) {
    overlap  <- as.integer(overlap)
    window.offset <- as.integer(window - round(window * (overlap /
        (overlap + 1))))
    chr.ordered <- unique(chromosomes)
    data.splitByChr <- split(data.frame(pos = positions, Bf, Af, good.reads),
        f = factor(chromosomes, levels = chr.ordered))
    lapply(data.splitByChr, function(data.oneChr) {
        range.pos <- range(data.oneChr$pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2] ) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        pos.cut <- cut(data.oneChr$pos, breaks = beam.coords)
        A.split <- split(data.oneChr$Af, f = pos.cut)
        B.split <- split(data.oneChr$Bf, f = pos.cut)
        d.split <- split(data.oneChr$good.reads, f = pos.cut)
        window.starts <- beam.coords[1:nWindows]
        window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        A.window <- lapply(idx.list, function(idx) {
            unlist(A.split[idx], use.names = FALSE)
        })
        B.window <- lapply(idx.list, function(idx) {
            unlist(B.split[idx], use.names = FALSE)
        })
        d.window <- lapply(idx.list, function(idx) {
            unlist(d.split[idx], use.names = FALSE)
        })
        window.quantiles <- mapply(b_allele_freq, Af = A.window, Bf = B.window,
            good.reads = d.window, conf = conf)
        window.counts <- sapply(B.window, length)
        data.frame(start = window.starts, end = window.ends,
            mean = window.quantiles[2, ], q0 = window.quantiles[1, ],
            q1 = window.quantiles[3, ], N = window.counts)
    })
}
