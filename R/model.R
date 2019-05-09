theoretical.depth.ratio <- function(CNt, cellularity, ploidy, CNn = 2,
    normal.ploidy = 2, avg.depth.ratio = 1) {
    cellu_copy_term <- (1 - cellularity) + (CNt / CNn * cellularity)
    ploidy_cellu_term <- (ploidy / normal.ploidy * cellularity) +
        1 - cellularity
    avg.depth.ratio * cellu_copy_term / ploidy_cellu_term
}

theoretical.baf <- function(CNt, B, cellularity, CNn = 2) {
    baf <- ( (B * cellularity) + ( 1 - cellularity) ) /
        ( (CNt * cellularity) + CNn * ( 1 - cellularity) )
    baf[CNn <= 1] <- NA
    baf
}

theoretical.mufreq <- function(CNt, Mt, cellularity, CNn = 2) {
    normal_alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
    all_alleles <- (CNt * cellularity) + CNn * (1 - cellularity)
    1 - (normal_alleles / all_alleles)
}

baf.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
    cn_ratio_vect <- seq(from = CNt.min / CNn, to = CNt.max / CNn,
        by = 1 / CNn)
    CNt <- cn_ratio_vect * CNn
    if (CNn < 2) {
        b_comb <- lapply(CNt, FUN = function(x) 0)
    } else {
        b_comb <- lapply(CNt, FUN = function(x) {
            seq(from = 0, to = trunc(x / 2))
        })
    }
    times_b <- sapply(b_comb, length)
    CNt <- rep(CNt, times = times_b)
    B <- unlist(b_comb)
    data.frame(CNn = CNn, CNt = CNt, B = B)
}

mufreq.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
    cn_ratio_vect <- seq(from = CNt.min / CNn,
        to = CNt.max / CNn, by = 1 / CNn)
    CNt <- cn_ratio_vect * CNn
    mut_comb <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
    times_muts <- sapply(mut_comb, length)
    data.frame(CNn = CNn, CNt = rep(CNt, times = times_muts),
        Mt = unlist(mut_comb))
}

baf.model.points <- function(cellularity, ploidy, baf_types, avg.depth.ratio) {
    depth_ratio <- theoretical.depth.ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = baf_types[, "CNn"], CNt = baf_types[, "CNt"],
        avg.depth.ratio = avg.depth.ratio)
    baf <- theoretical.baf(cellularity = cellularity, CNn = baf_types[, "CNn"],
        CNt = baf_types[, "CNt"], B = baf_types[, "B"])
    data.frame(BAF = baf, depth.ratio = depth_ratio)
}

mufreq.model.points <- function(cellularity, ploidy, mufreq_types,
    avg.depth.ratio) {
    mufreqs <- theoretical.mufreq(cellularity = cellularity,
        CNn = mufreq_types[, "CNn"], CNt = mufreq_types[, "CNt"],
        Mt = mufreq_types[, "Mt"])
    depth_ratio <- theoretical.depth.ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = mufreq_types[, "CNn"],
        CNt = mufreq_types[, "CNt"], avg.depth.ratio = avg.depth.ratio)
    data.frame(mufreqs = mufreqs, depth.ratio = depth_ratio)
}

b_allele_freq <- function(Af, Bf, good.reads, conf = 0.95) {
    if (length(Bf) > 1) {
        dd <- density(c(Bf, Af), weight = c(good.reads, good.reads) /
            (2 * sum(good.reads)))
        points.max <- which(diff(sign(diff(dd$y))) == -2) + 1
        if (length(points.max) < 1) {
            points.max <- which(diff(sign(diff(dd$y))) == -1) + 1
        }
        l.max <- dd$x[points.max]
        d.max <- dd$y[points.max]
        b.val <- l.max[which.max(dd$y[dd$x %in% l.max])]
        if (length(b.val) < 1) {
            message('WARNING', l.max, d.max, b.val)
            b.val <- min(l.max)
        }
        d.val <- d.max[which(l.max == b.val)]
        b.range <- range(dd$x[dd$y >=  d.val - (d.val * (1 - conf))])
        if (b.val > 0.5) {
            b.val <- 1 - b.val
        }
        max_diff <- max(b.range) - b.val
        min_diff <- b.val - min(b.range)
        min_diff <- min(c(max_diff, min_diff))
        c(b.val - min_diff, b.val, b.val + min_diff)
    } else if (length(Bf) == 1) {
        c(Bf, Bf, Bf)
    } else {
        c(NA, NA, NA)
    }
}
