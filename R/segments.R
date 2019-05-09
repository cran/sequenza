extract_breaks <- function(data, data_het, ratio, baf, breaks, gamma, kmin,
    gamma.pcf, kmin.pcf, assembly, chromosome,
    method = c("het", "full", "fast")) {
    method_list <- c("het", "full", "fast")
    if (is.null(breaks)) {
        if (method %in% method_list) {
            if (method == "fast"){
                breaks <- breaks_fast(ratio_win = ratio, baf_win = baf,
                    gamma = gamma, kmin = kmin, chr = chromosome)
            } else {
                breaks <- breaks_het(data = data_het, gamma = gamma,
                    kmin = kmin, assembly = assembly)
            }
            if (method == "full") {
                breaks <- breaks_full(data = data, gamma = gamma.pcf,
                    kmin = kmin.pcf, assembly, breaks.het = breaks)
            }
        } else {
            stop("Available methods are \'full\', \'het\' and \'fast\'.")
        }
    }
    breaks
}


breaks_het <- function(data, gamma, kmin, assembly){
    try(
        find.breaks(data, gamma = gamma, assembly = assembly,
            kmin = kmin, baf.thres = c(0, 0.5)),
    silent = FALSE)
}

breaks_full <- function(data, gamma, kmin,
    assembly, breaks.het = NULL) {
    merge.breaks <- function (breaks, breaks.het) {
        merged.breaks <- unique(sort(c(breaks$start.pos,
            breaks$end.pos, breaks.het$start.pos, breaks.het$end.pos)))
        merged.breaks <- merged.breaks[diff(merged.breaks) > 1]
        merged.start <- merged.breaks
        merged.start[-1] <- merged.start[-1] + 1
        breaks <- data.frame(chrom = unique(breaks$chrom),
            start.pos = merged.start[- (length(merged.start))],
            end.pos = merged.breaks[-1])
    }
    breaks <- find.breaks(data, gamma = gamma, kmin = kmin,
        assembly = assembly, seg.algo = "pcf")
    if (!is.null(breaks.het)) {
        chr.p <- merge.breaks(breaks[breaks$arm == "p", ],
            breaks.het[breaks.het$arm == "p", ])
        chr.q <- merge.breaks(breaks[breaks$arm == "q", ],
            breaks.het[breaks.het$arm == "q", ])
        breaks <- rbind(chr.p, chr.q)
    }
    breaks
}

breaks_fast <- function(ratio_win, baf_win, chr, gamma, kmin) {
    BAF <- data.frame(chrom = chr,
        pos = c(baf_win[[1]]$start, tail(baf_win[[1]]$end, n = 1)),
        s1 = c(baf_win[[1]]$mean, tail(baf_win[[1]]$mean, n = 1)))
    logR <- data.frame(chrom = chr,
        pos = c(ratio_win[[1]]$start, tail(ratio_win[[1]]$end, n = 1)),
        s1 = c(log2(ratio_win[[1]]$mean),
            log2(tail(ratio_win[[1]]$mean, n = 1))))
    cat(nrow(BAF), nrow(logR), "\n")
    not.cover <- is.na(logR$s1) & is.na(BAF$s1)
    BAF  <- BAF[!not.cover, ]
    logR <- logR[!not.cover, ]
    logR.wins <- copynumber::winsorize(logR, verbose = FALSE)
    allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF,
        baf.thres = c(0, 0.5), verbose = FALSE,
        gamma = gamma, kmin = kmin)
    if (length(grep("^chr", chr)) > 0) {
        allele.seg$chrom <- paste0("chr", allele.seg$chrom)
    }
    breaks   <- allele.seg[, c("chrom", "start.pos", "end.pos")]
    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 0))
    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
    breaks
}


find.breaks <- function(seqz.baf, gamma = 80, kmin = 10,
    baf.thres = c(0, 0.5), verbose = FALSE, seg.algo = "aspcf", ...) {
    chromosome <- gsub(x = seqz.baf$chromosome,
        pattern = "chr", replacement = "")
    logR = data.frame(chrom = chromosome,
        pos = seqz.baf$position,
        s1 = log2(seqz.baf$adjusted.ratio))
    logR.wins <- copynumber::winsorize(logR, verbose = verbose)
    if (seg.algo == "aspcf"){
        BAF = data.frame(chrom = chromosome,
            pos = seqz.baf$position,
            s1 = seqz.baf$Bf)
        allele.seg <- copynumber::aspcf(logR = logR.wins,
            BAF = BAF, baf.thres = baf.thres, verbose = verbose,
            gamma = gamma, kmin = kmin, ...)
    } else if (seg.algo == "pcf") {
        allele.seg <- copynumber::pcf(data = logR.wins, verbose = verbose,
            gamma = gamma, kmin = kmin, ...)
    } else {
      stop("Segmentation algorithm must be either \'aspcf\' or \'pcf\'.")
    }
    if (length(grep("chr", seqz.baf$chromosome)) > 0) {
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
    breaks   <- allele.seg[, c("chrom", "start.pos", "end.pos", "arm")]
    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1],0))
    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
    breaks
}

segment.breaks <- function(seqz.tab, breaks, min.reads.baf = 1,
    weighted.mean = TRUE) {
    if (weighted.mean){
        w.r <- sqrt(seqz.tab$depth.normal)
        rw <- seqz.tab$adjusted.ratio * w.r
        w.b <- sqrt(seqz.tab$good.reads)
        bw <- seqz.tab$Bf * w.b
        seqz.tab <- cbind(seqz.tab[, c("chromosome", "position",
            "zygosity.normal", "good.reads", "Af", "Bf")],
            rw = rw, w.r = w.r, bw = bw, w.b = w.b)
    }
    chr.order <- unique(seqz.tab$chromosome)
    seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
    segments <- list()
    for (i in 1:length(seqz.tab)) {
        seqz.b.i <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == "het", ]
        seqz.b.i <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, ]
        breaks.i <- breaks[breaks$chrom == names(seqz.tab)[i], ]
        nb <- nrow(breaks.i)
        breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,
            c("start.pos", "end.pos")], f = 1:nb))
        unique.breaks <- function(b, offset = 1) {
            while(any(diff(b) == 0)) {
                b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset
            }
            b
        }
        breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)
        fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
        fact.b.i <- cut(seqz.b.i$position, breaks.vect)
        seg.i.s.r <- sapply(X = split(seqz.tab[[i]]$chromosome,
            f = fact.r.i), FUN = length)
        seg.i.s.b <- sapply(X = split(seqz.b.i$chromosome,
            f = fact.b.i), FUN = length)

        if (weighted.mean) {
            seg.i.rw    <- sapply(X = split(seqz.tab[[i]]$rw, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.w.r   <- sapply(X = split(seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.r.sd  <- sapply(X = split(seqz.tab[[i]]$rw /
                seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd  <- sapply(X = split(seqz.b.i$bw /
                seqz.b.i$w.b, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)
            window.quantiles <- mapply(b_allele_freq, Af = A.split,
                Bf = B.split, good.reads = d.split, conf = 0.95)
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw / seg.i.w.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        } else {
            seg.i.r    <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) mean(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)
            window.quantiles <- mapply(b_allele_freq, Af = A.split,
                Bf = B.split, good.reads = d.split, conf = 0.95)
            seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        }
        segments[[i]] <- segments.i[seq(from = 1,
            to = nrow(segments.i), by = 2),]
    }
    segments <- do.call(rbind, segments[as.factor(chr.order)])
    row.names(segments) <- 1:nrow(segments)
    len.seg <- (segments$end.pos - segments$start.pos) / 1e6
    segments[(segments$N.ratio / len.seg) >= 2, ]
}
