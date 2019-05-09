plotWindows <- function(seqz.window, m.lty = 1, m.lwd = 3,
                         m.col = "black", q.bg = "lightblue",
                         log2.plot = FALSE, n.min = 1, xlim, ylim,
                         add = FALSE, ...) {
    if (log2.plot) {
        seqz.window[, c(3, 4, 5)] <- log2(seqz.window[, c(3, 4, 5)])
    }
    if (!add) {
        if (missing(xlim))
            xlim <- c(seqz.window$start[1], seqz.window$end[nrow(seqz.window)])
        if (missing(ylim))
           ylim <- c(min(seqz.window$q0, na.rm = TRUE),
               max(seqz.window$q1, na.rm = TRUE))
        plot(xlim, ylim, type = "n", ...)
    }
    seqz.window <- seqz.window[seqz.window$N >= n.min, ]
    rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
        xright = seqz.window$end, ytop = seqz.window$q1,
        col = q.bg, border = NA)
    segments(y0 = seqz.window$mean, x0 = seqz.window$start,
        x1 = seqz.window$end, lty = m.lty, lwd = m.lwd, col = m.col)

}

gc.plot <- function(gc_list, range.gc = NULL, range.depth = NULL, ...) {
    n <- gc_list$n
    n[n == 0] <- NA
    gc <- gc_list$gc
    depth <- gc_list$depth
    if (length(range.gc) == 2) {
        range.gc <- sort(range.gc)
        gc.select <- gc <= range.gc[2] & gc >= range.gc[1]
        gc <- gc[gc.select]
        n <- n[gc.select, ]
    }
    if (length(range.depth) == 2) {
        range.depth <- sort(range.depth)
        depth.select <- depth <= range.depth[2] & depth >= range.depth[1]
        depth <- depth[depth.select]
        n <- n[, depth.select]
    }
    colorgram(x = gc, y = depth, z = n,  ...)
}

gc.summary.plot <- function(gc_list, mean.col = 1, median.col = 2,
    scale.subset = 1.5, ...){
    mengc <- mean_gc(gc_list)
    medgc <- median_gc(gc_list)
    max_depth <- round(max(c(mengc, medgc)) * scale.subset, 0)
    gc.plot(gc_list, range.depth = c(0, max_depth), ...)
    lines(x = gc_list$gc, y = mengc, lwd = 3, col = mean.col)
    lines(x = gc_list$gc, y = medgc, lwd = 3, col = median.col)
    legend("topright", c("Mean depth", "Median depth"),
        col = c(mean.col, median.col), bg = "white", lty = 1, lwd = 3)
}

cp.plot <- function (cp.table, xlab = "Ploidy", ylab = "Cellularity",
    zlab = "Scaled rank LPP",
    colFn = colorRampPalette(c("white", "lightblue")), ...) {
    z <- matrix(rank(cp.table$lpp), nrow = nrow(cp.table$lpp)) /
        length(cp.table$lpp)
    map <- makecmap(c(0, 1), colFn = colFn, include.lowest = TRUE)
    colorgram(x = cp.table$ploidy, y = cp.table$cellularity, z = z,
        map = map, las = 1, xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

cp.plot.contours <- function(cp.table, likThresh = c(0.95),
    alternative = TRUE, col = palette(), legend.pos = "bottomright",
    pch = 18, alt.pch = 3, ...) {
    znormsort <- sort(cp.table$lpp, decreasing = TRUE)
    znormcumLik <- cumsum(znormsort)
    n <- sapply(likThresh, function(x) sum(znormcumLik < x) + 1)
    LikThresh <- znormsort[n]
    names(LikThresh) <- paste0(likThresh * 100, "%")
    contour(x = cp.table$ploidy, y = cp.table$cellularity, z = cp.table$lpp,
        levels = znormsort[n], col = col, drawlabels = FALSE,
        xlab = "Ploidy", ylab = "Cellularity", ...)
    max.xy <- which(cp.table$lpp == max(cp.table$lpp), arr.ind = TRUE)
    points(x = cp.table$ploidy[max.xy[, "row"]],
        y = cp.table$cellularity[max.xy[, "col"]], pch = pch)
    if (alternative == TRUE){
        alt.sol <- alternative.cp.solutions(cp.table)
        alt.sol <- alt.sol[-1, ]
        points(x = alt.sol$ploidy, y = alt.sol$cellularity, pch = alt.pch)
    }
    if (!is.na(legend.pos)) {
        if (alternative == FALSE) {
            legend(legend.pos, legend = c(paste("C.R.", names(LikThresh),
                sep = " "), "Point estimate"),
                col = c(col[1:length(LikThresh)], "black"),
                lty = c(rep(1, length(LikThresh)), NA),
                pch = c(rep(NA, length(LikThresh)), pch),
                border = NA, bty = "n")
        } else {
            legend(legend.pos, legend = c(paste("C.R.",
                    names(LikThresh), sep = " "),
                    "Point estimate", "Alternative solutions"),
                col = c(col[1:length(LikThresh)], "black", "black"),
                lty = c(rep(1, length(LikThresh)), NA, NA),
                pch = c(rep(NA, length(LikThresh)), pch, alt.pch),
                border = NA, bty = "n")
        }
    }
   invisible(LikThresh)
}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL,
    segments = NULL,  min.N.baf = 1, min.N.ratio = 1e4, main = "",
    vlines = FALSE, legend.inset = c(-20 * strwidth("a", units = "figure"), 0),
    CNn = 2, cellularity = NULL, ploidy = NULL,
    avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24", model.col = 1,
    x.chr.space = 10) {
    if (is.null(segments)) {
      data.model <- NULL
    } else {
        if ("CNt" %in% colnames(segments)) {
            if (length(c(cellularity, ploidy, avg.depth.ratio)) != 3) {
                data.model <- NULL
            } else {
                data.model <- list()
                CNt.max <- max(segments$CNt, na.rm = TRUE) + 1
                CNt.min <- 0
                baf_types <- baf.types.matrix(CNt.min = CNt.min,
                    CNt.max = CNt.max, CNn = 2)
                data.model$baf <- baf.model.points(cellularity = cellularity,
                    ploidy = ploidy, baf_types = baf_types,
                    avg.depth.ratio = avg.depth.ratio)
                data.model$baf <- data.frame(CNt =  baf_types$CNt,
                    A = baf_types$CNt - baf_types$B, B = baf_types$B,
                    data.model$baf)
                mufreq_types <- mufreq.types.matrix(CNt.min = CNt.min,
                    CNt.max = CNt.max, CNn = CNn)
                data.model$muf <- cbind(mufreq_types,
                    mufreq.model.points(cellularity = cellularity,
                        ploidy = ploidy, mufreq_types = mufreq_types,
                        avg.depth.ratio = avg.depth.ratio))
            }
        } else {
            data.model <- NULL
        }
    }
    if (is.null(mut.tab)) {
        par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),
            mfcol = c(2, 1), xaxt = "n")
        min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
        max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
        xlim <- c(min.x, max.x)
    } else {
        min.x <- min(c(min(baf.windows$start), min(ratio.windows$start),
            min(mut.tab$position)))
        max.x <- max(c(max(baf.windows$end), max(ratio.windows$end),
            max(mut.tab$position)))
        xlim <- c(min.x, max.x)
        par(mar = c(0, 4, 0, 10), oma = c(5, 0, 4, 0),
            mfcol = c(3, 1), xaxt = "n", xpd = TRUE)
        mutation.colors <- c(
            "A>C" = rgb(red = 0, green = 178, blue = 238,
                alpha = 120, maxColorValue = 255),
            "T>G" = rgb(red =   0, green = 178, blue = 238,
                alpha = 120, maxColorValue = 255),
            "A>G" = rgb(red = 255, green = 64, blue = 64,
                alpha = 120, maxColorValue = 255),
            "T>C" = rgb(red = 255, green = 64, blue = 64,
                alpha = 120, maxColorValue = 255),
            "A>T" = rgb(red =  34, green = 139, blue = 34,
                alpha = 120, maxColorValue = 255),
            "T>A" = rgb(red =  34, green = 139, blue = 34,
                alpha = 120, maxColorValue = 255),
            "C>A" = rgb(red = 139, green = 90, blue = 0,
                alpha = 120, maxColorValue = 255),
            "G>T" = rgb(red = 139, green = 90, blue = 0,
                alpha = 120, maxColorValue = 255),
            "C>G" = rgb(red = 127, green =   0, blue = 255,
                alpha = 120, maxColorValue = 255),
            "G>C" = rgb(red = 127, green =   0, blue = 255,
                alpha = 120, maxColorValue = 255),
            "C>T" = rgb(red = 255, green = 215, blue = 0,
                alpha = 120, maxColorValue = 255),
            "G>A" = rgb(red = 255, green = 215, blue = 0,
                alpha = 120, maxColorValue = 255))
        plot(x = mut.tab$position, y = mut.tab$F,
            ylab = "Mutant allele frequency", las = 1, pch = 19,
            col = c(mutation.colors, "NA" = NA)[as.character(mut.tab$mutation)],
            ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
        unique.colors <- unique(mutation.colors)
        labels <- sapply(unique.colors, function(a) {
            paste(names(mutation.colors)[mutation.colors == a],
                collapse = ", ")
            })
        legend(y = "center", x  = "right", legend = labels,
            inset = legend.inset, pch = 19, col = unique.colors,
            pt.bg = unique.colors, border = NA, bty = "n")
        if (!is.null(segments)){
            if (vlines) {
                abline(v = segments$end.pos, lwd = 1, lty = 2)
            }
            if (!is.null(data.model)) {
                for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i],
                    x1 = segments$end.pos[i],
                    y0 = unique(data.model$muf$mufreqs[
                        data.model$muf$CNt == segments$CNt[i]]),
                    lwd = model.lwd, lty = model.lty, col = model.col)
                }
            }
        }
    }
    if (!is.null(segments)){
        plot(ylab = "B allele frequency", type = "n",
            x = xlim, y = c(0, 0.5), las = 1)
        plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf, add = TRUE)
        if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
        }
        segments(x0 = segments$start.pos, y0 = segments$Bf,
            x1 = segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
        if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
               segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                        y0 = unique(data.model$baf$BAF[
                            data.model$baf$CNt == segments$CNt[i]]),
                        lwd = model.lwd, lty = model.lty, col = model.col)
            }
        }
    } else {
        plotWindows(baf.windows, ylab = "B allele frequency",
            xlim = xlim, ylim = c(0, 0.5), las = 1, n.min = min.N.baf)
    }
    plotWindows(ratio.windows, ylab = "Depth ratio",
        las = 1, n.min = min.N.ratio, ylim = c(0, 2.5))
    if (!is.null(segments)){
        if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
        }
        segments(x0 = segments$start.pos, y0 = segments$depth.ratio,
            x1 = segments$end.pos, y1 = segments$depth.ratio,
            col = "red", lwd = 3)
        if (!is.null(data.model)) {
            ratios.theoric <- unique(data.model$muf[, c("CNt", "depth.ratio")])
            segments(x0 = rep(min(segments$start.pos, na.rm = TRUE),
                    times = nrow(ratios.theoric)),
                x1 = rep(max(segments$end.pos, na.rm = TRUE),
                    times = nrow(ratios.theoric)),
                y0 = ratios.theoric$depth.ratio, lwd = model.lwd,
                lty = model.lty, col = model.col)

            axis(labels = as.character(ratios.theoric$CNt), side = 4,
                line = 0, las = 1, at = ratios.theoric$depth.ratio)
            mtext(text = "Copy number", side = 4, line = 2,
                cex = par("cex.lab") * par("cex"))
        }
    }
    par(xaxt = "s")
    axis(labels = as.character(round(seq(xlim[1] / 1e6, xlim[2] / 1e6,
            by = x.chr.space), 0)),
        side = 1, line = 0, at = seq(xlim[1], xlim[2], by = 1e6 * x.chr.space),
        outer = FALSE, cex = par("cex.axis") * par("cex"))
    mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
        cex = par("cex.lab") * par("cex"))
    mtext(main, 3, outer = TRUE, cex = par("cex.main") * par("cex"), line = 2)
}

genome.view <- function(seg.cn, info.type = "AB", ...) {
    chr.order <- unique(seg.cn$chromosome)
    seg.list <- split(x = seg.cn[,
        c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")],
        f = seg.cn$chromosome)

    seg.list <- seg.list[order(order(chr.order))]
    seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
    seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
    seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
    chr.offset <- 0
    for (i in 1:length(seg.pos)){
        seg.pos[[i]] <- seg.pos[[i]] + chr.offset
        colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
        chr.offset <- seg.max[i]
    }
    seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
    abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
    abs.segments <- do.call(rbind, abs.list)
    if (info.type == "AB") {
        na_As <- is.na(abs.segments$A)
        max_A <- max(abs.segments$A, na.rm = TRUE)
        abs.segments$A[na_As] <- abs.segments$CNt[na_As]
        plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
            y = c(-0.1, (max_A + 0.1)), type = "n",
            ylab = "Copy number", xlab = "Position (Mb)",
            xaxt = "n",  yaxt = "n", xaxs = "i", ...)
        axis(labels = 0:max_A, at = 0:max_A, side = 2, line = 0, las = 1)
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = (abs.segments$B - 0.1), y1 = (abs.segments$B - 0.1),
            col = "blue", lwd = 5, lend = 1)
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = (abs.segments$A + 0.1), y1 = (abs.segments$A + 0.1),
            col = "red", lwd = 5, lend = 1)
    } else {
        min_CNt <- min(abs.segments$CNt, na.rm = TRUE)
        max_CNt <- max(abs.segments$CNt, na.rm = TRUE)
        plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
            y = c(min_CNt, max_CNt), type = "n",
            ylab = "Copy number", xlab = "Position (Mb)",
            xaxt = "n", yaxt = "n", xaxs = "i", ...)
        axis(labels = min_CNt:max_CNt,
            at = min_CNt:max_CNt,
            side = 2, line = 0, las = 1)
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = abs.segments$CNt, y1 = abs.segments$CNt, col = "red",
            lwd = 5, lend = 1)
    }
    abline(v = c(0, seg.max), lty = 3)
    for (i in 1:length(abs.list)){
        max.pos <- nrow(abs.list[[i]])
        mtext(chr.order[i], side = 3, line = 0,
            at = sum(abs.list[[i]]$abs.start[1],
                abs.list[[i]]$abs.end[max.pos]) / 2)
    }
    axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / 1e6,
        abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / 1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1],
            abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7),
        outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
        line = 1)
}

plotRawGenome <- function(sequenza.extract, cellularity,
    ploidy, CNt.max = 7, main = "", mirror.BAF = TRUE, ...){
    max.end <- sapply(sequenza.extract$ratio, FUN = function(x) {
        max(x$end, na.rm = TRUE)
    })
    max.end <- c(0, cumsum(as.numeric(max.end)))
    chrs <- names(sequenza.extract$ratio)
    coords.names <- (max.end + c(diff(max.end) / 2, 0))[1:length(chrs)]
    new.coords <- function(win.list, max.end){
        lapply(1:length(win.list), FUN = function(x) {
            y <- win.list[[x]]
            y$start <- y$start + max.end[x]
            y$end <- y$end + max.end[x]
            y
        }
    )}
    new.coords.segs <- function(segs, max.end){
        lapply(1:length(segs), FUN = function(x) {
            y <- segs[[x]]
            y$start.pos <- y$start.pos + max.end[x]
            y$end.pos <- y$end.pos + max.end[x]
            y
        }
    )}

    ratio.new <- new.coords(sequenza.extract$ratio, max.end)
    BAF.new <- new.coords(sequenza.extract$BAF, max.end)

    segs.new <- do.call(rbind,
        new.coords.segs(sequenza.extract$segments, max.end))
    avg.depth.ratio <- 1

    par(mar = c(1, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2, 1), ...)

    if (mirror.BAF) {
        AAF.new <- lapply(BAF.new, function(x) {
            x[, 3:5] <- 1 - x[, 3:5]
            x
        })
        plot(x = c(min(max.end), max(max.end)), y = c(0, 1), main = main,
            xlab = NA, ylab = "Allele frequency", type = "n", las = 1,
            xaxs = "i", yaxs = "i", xaxt = "n" )
            plotWindows(seqz.window = do.call(rbind, AAF.new),
                q.bg = "lightblue", m.col = "black", add = T)
            segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
                y0 = 1 - (segs.new$Bf), y1 = 1 - (segs.new$Bf),
                col = "red", lwd = 2, lend = 1)
    } else {
        plot(x = c(min(max.end), max(max.end)), y = c(0, 0.5), main = main,
            xlab = NA, ylab = "B allele frequency", type = "n", las = 1,
            xaxs = "i", yaxs = "i", xaxt = "n" )
    }
    plotWindows(seqz.window = do.call(rbind, BAF.new), q.bg = "lightblue",
        m.col = "black", add = T)
    segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
        y0 = segs.new$Bf, y1 = segs.new$Bf, col = "red", lwd = 2, lend = 1)
    abline(v = max.end, lty = 1)
    plot(x = c(min(max.end), max(max.end)), y = c(0, 2.5), main = "",
        xlab = NA, ylab = "Depth ratio", type = "n", las = 1, xaxs = "i",
        yaxs = "i", xaxt = "n")
    plotWindows(seqz.window = do.call(rbind, ratio.new), q.bg = "lightblue",
        m.col = "black", add = T)
    segments(x0 = segs.new$start.pos, x1 = segs.new$end.pos,
        y0 = (segs.new$depth.ratio), y1 = (segs.new$depth.ratio),
        col = "red", lwd = 2, lend = 1)
    if (!missing(ploidy) & !missing(cellularity)){
        types <- baf.types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2)
        depth.ratios <- baf.model.points(cellularity = cellularity,
            ploidy = ploidy, avg.depth.ratio = avg.depth.ratio,
            baf_types = types)[, "depth.ratio"]
        depth.ratios <- unique(data.frame(CNt = types$CNt,
            ratio = depth.ratios))
        abline(h = depth.ratios$ratio, lty = 2)
        axis(labels = as.character(depth.ratios$CNt), side = 4,
            line = 0, las = 1, at = depth.ratios$ratio)
        mtext(text = "Copy number", side = 4, line = 2,
            cex = par("cex.lab") * par("cex"))
    }
    abline(v = max.end, lty = 1)
    axis(labels = chrs, at = coords.names, side = 1, cex.axis = 1)
}

baf.model.view <- function(cellularity, ploidy, segs,
    BAF.space = seq(0.001, 0.5, 0.005), ratio.space = seq(0.01, 2.5, 0.05),
    avg.depth.ratio = 1, CNt.max = 7, segment.filter = 3e6, col = "black") {
    s.b   <- mean(segs$sd.BAF, na.rm = TRUE)
    s.r   <- mean(segs$sd.ratio, na.rm = TRUE)
    l.s   <- segs$end.pos - segs$start.pos
    s.big <- l.s >= segment.filter
    test.values <- expand.grid(Bf = BAF.space, ratio = ratio.space,
        KEEP.OUT.ATTRS = FALSE)
    both.space  <- baf.bayes(Bf = test.values$Bf, CNt.max = CNt.max,
        CNt.min = 0, depth.ratio = test.values$ratio,
        cellularity = cellularity, ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio, sd.Bf = s.b, weight.Bf = 10,
        sd.ratio = s.r, weight.ratio = 10, ratio.priority = F, CNn = 2)
    both.space <- as.data.frame(both.space)
    z <- tapply(both.space$LPP, list(test.values$Bf, test.values$ratio), mean)
    x <- as.numeric(rownames(z))
    y <- as.numeric(colnames(z))
    t <- baf.types.matrix(CNt.min = 0, CNt.max = CNt.max, CNn = 2)
    mpts <- cbind(t, baf.model.points(cellularity = cellularity,
        ploidy = ploidy, baf_types = t, avg.depth.ratio = avg.depth.ratio)
    )
    mpts <- unique(mpts[, c("CNt", "depth.ratio")])
    par(mar = c(5.1, 4.1, 4.1, 4.1))
    rev.heat <- function(...){rev(heat.colors(...))}
    suppressWarnings(colorgram(x, y, z, key = NA, nz = 1000,
        xlab = "B allele frequency", ylab = "Depth ratio",
        main = paste("cellularity:", cellularity, "ploidy:", ploidy,
        "sd.BAF:", round(s.b, 2), sep = " "),
        map = makecmap(z, breaks = unique(quantile(z, seq(0.25, 1, 0.0001))),
            right = TRUE, n = 1000, colFn = rev.heat), outlier = "white",
        las = 1, xlim = c(0, 0.5)))
    axis(side = 4, at = mpts$depth.ratio, labels = mpts$CNt, las = 1)
    mtext(text = "Copy number", side = 4, line = 2)
    segs$col <- col
    points(x = segs$Bf[s.big], y = segs$depth.ratio[s.big],
        pch = 1, cex = 1, col = segs$col[s.big])
    points(x = segs$Bf[!s.big], y = segs$depth.ratio[!s.big], pch = ".",
        cex = 1, col = segs$col[!s.big])
}
