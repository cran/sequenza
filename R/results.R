sequenza.results <- function(sequenza.extract, cp.table = NULL,
    sample.id, out.dir = getwd(), cellularity = NULL, ploidy = NULL,
    female = TRUE, CNt.max = 20, ratio.priority = FALSE,
    XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
    if(!file.exists(out.dir)) {
        dir.ok <- dir.create(path = out.dir, recursive = TRUE)
        if(!dir.ok) {
            stop("Directory does not exist and cannot be created: ", out.dir)
        }
    }
    make_filename <- function(x){
        file.path(out.dir, paste(sample.id, x, sep = "_"))
    }
    cp.file   <- make_filename("CP_contours.pdf")
    cint.file <- make_filename("confints_CP.txt")
    chrw.file <- make_filename("chromosome_view.pdf")
    depths.file <- make_filename("chromosome_depths.pdf")
    gc.file <- make_filename("gc_plots.pdf")
    geno.file <- make_filename("genome_view.pdf")
    cn.file <- make_filename("CN_bars.pdf")
    fit.file <- make_filename("model_fit.pdf")
    alt.file <- make_filename("alternative_solutions.txt")
    afit.file <- make_filename("alternative_fit.pdf")
    muts.file <- make_filename("mutations.txt")
    segs.file <- make_filename("segments.txt")
    robj.extr <- make_filename("sequenza_extract.RData")
    robj.fit <- make_filename("sequenza_cp_table.RData")
    log.file <- make_filename("sequenza_log.txt")
    seg.tab <- do.call(rbind, sequenza.extract$segments[chromosome.list])
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos) / 1e6

    avg.depth.ratio <- sequenza.extract$avg.depth.ratio
    assign(x = paste0(sample.id, "_sequenza_extract"),
        value = sequenza.extract)
    save(list = paste0(sample.id, "_sequenza_extract"), file = robj.extr)
    if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))){
        stop("cp.table and/or cellularity and ploidy argument are required.")
    }

    pdf(gc.file, width = 10, height = 5)
    par(mfrow=c(1, 2))
    gc.summary.plot(sequenza.extract$gc$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc_norm$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the tumor sample")
    gc.summary.plot(sequenza.extract$gc_norm$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the tumor sample")
    dev.off()
    pdf(depths.file, height = 10, width = 15)
    for (i in unique(seg.tab$chromosome)) {
        max_coord_chr_i <- max(sequenza.extract$ratio[[i]]$end)
        par(mfcol = c(3, 2), xaxt = "n", mar = c(0, 4, 3, 0), oma = c(5, 0, 4, 0))
        plotWindows(sequenza.extract$depths$raw$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5),
            main = paste("raw", i, sep = " "))
        plotWindows(sequenza.extract$depths$raw$tumor[[i]],
            ylab = "tumor depth", ylim = c(0, 2.5))
        plotWindows(sequenza.extract$raw_ratio[[i]],
            ylab = "depth ratio", ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i / 1e6,
                by = 10), 0)),
            side = 1, line = 0, at = seq(0, max_coord_chr_i, by = 1e7),
            outer = FALSE, cex = par("cex.axis") * par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
        par(xaxt = "n")
        plotWindows(sequenza.extract$depths$norm$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5),
            main = paste("normalized", i, sep = " "))
        plotWindows(sequenza.extract$depths$norm$tumor[[i]],
            ylab = "tumor depth", ylim = c(0, 2.5))
        plotWindows(sequenza.extract$ratio[[i]],
            ylab = "depth ratio", ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i / 1e6,
                by = 10), 0)),
            side = 1, line = 0, at = seq(0, max_coord_chr_i, by = 1e7),
            outer = FALSE, cex = par("cex.axis") * par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
    }
    dev.off()

    if (!is.null(cp.table)){
        assign(x = paste0(sample.id, "_sequenza_cp_table"), value = cp.table)
        save(list = paste0(sample.id, "_sequenza_cp_table"), file = robj.fit)
        cint <- get.ci(cp.table)
        pdf(cp.file)
            cp.plot(cp.table)
            cp.plot.contours(cp.table, add = TRUE,
                likThresh = c(0.95), col = "red", pch = 20)
            if (!is.null(cellularity) || !is.null(ploidy)) {
                if (is.null(cellularity)) {
                    cellularity <- cint$max.cellularity
                }
                if (is.null(ploidy)) {
                    ploidy <- cint$max.ploidy
                }
                points(x = ploidy, y = cellularity, pch = 5)
                text(x = ploidy, y = cellularity, labels = "User selection",
                    pos = 3, offset = 0.5)
            } else {
                cellularity <- cint$max.cellularity
                ploidy <- cint$max.ploidy
            }
        dev.off()
    }
    mut.tab <- na.exclude(do.call(rbind,
        sequenza.extract$mutations[chromosome.list]))
    if (female){
        segs.is.xy <- seg.tab$chromosome == XY["Y"]
        mut.is.xy <- mut.tab$chromosome == XY["Y"]
    } else{
        segs.is.xy <- seg.tab$chromosome %in% XY
        mut.is.xy <- mut.tab$chromosome %in% XY
    }
    avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE) /
        sum(seg.tab$N.ratio, na.rm = TRUE)
    avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE) /
        sum(seg.tab$N.BAF, na.rm = TRUE)
    cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
        depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
        cellularity = cellularity, ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio,
        sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
        weight.ratio = seg.len[!segs.is.xy],
        sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
        weight.Bf = 1, ratio.priority = ratio.priority, CNn = 2)
    seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
    if (!female){
        if (sum(segs.is.xy) >= 1) {
            cn.alleles  <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio,
                sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                weight.ratio = seg.len[segs.is.xy], sd.Bf = NA,
                weight.Bf = NA, ratio.priority = ratio.priority, CNn = 1)
            seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
            seg.res <- rbind(seg.res, seg.xy)
        }
    }
    write.table(seg.res, file = segs.file, col.names = TRUE,
        row.names = FALSE, sep = "\t", quote = FALSE)
    if (nrow(mut.tab) > 0) {
        mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy],
            CNt.max = CNt.max,
            depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
            cellularity = cellularity, ploidy = ploidy,
            avg.depth.ratio = avg.depth.ratio, CNn = 2)
        mut.res <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
        if (!female){
            if (sum(mut.is.xy) >= 1) {
                mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy],
                    CNt.max = CNt.max,
                    depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                    cellularity = cellularity, ploidy = ploidy,
                    avg.depth.ratio = avg.depth.ratio, CNn = 1)
                mut.xy <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
                mut.res <- rbind(mut.res, mut.xy)
            }
        }
        write.table(mut.res, file = muts.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }
    pdf(chrw.file)
    for (i in unique(seg.res$chromosome)) {
        if (!female && i %in% XY){
            CNn <- 1
        } else {
            CNn <- 2
        }
        chromosome.view(mut.tab = sequenza.extract$mutations[[i]],
            baf.windows = sequenza.extract$BAF[[i]],
            ratio.windows = sequenza.extract$ratio[[i]],
            cellularity = cellularity, ploidy = ploidy, main = i,
            segments = seg.res[seg.res$chromosome == i, ],
            avg.depth.ratio = avg.depth.ratio, CNn = CNn, min.N.ratio = 1)
    }
    dev.off()
    pdf(geno.file, height = 5, width = 15)
    if (sum(!is.na(seg.res$A)) > 0) {
        genome.view(seg.res)
    }
    genome.view(seg.res, "CN")
    plotRawGenome(sequenza.extract, cellularity = cellularity, ploidy = ploidy,
        mirror.BAF = TRUE)
    dev.off()
    barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos,
        CNt = seg.res$CNt)
    cn.sizes <- split(barscn$size, barscn$CNt)
    cn.sizes <- sapply(cn.sizes, "sum")
    pdf(cn.file)
    barplot(round(cn.sizes / sum(cn.sizes) * 100), names = names(cn.sizes),
        las = 1, ylab = "Percentage (%)", xlab = "Copy number")
    dev.off()

    ## Write down the results.... ploidy etc...
    if (!is.null(cp.table)){
        res.tab <- data.frame(cellularity = c(cint$confint.cellularity[1],
            cint$max.cellularity[1], cint$confint.cellularity[2]),
            ploidy.estimate = c(cint$confint.ploidy[1], cint$max.ploidy[1],
                cint$confint.ploidy[2]),
            ploidy.mean.cn = weighted.mean(x = as.integer(names(cn.sizes)),
                w = cn.sizes))
        write.table(res.tab, cint.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }
    pdf(fit.file, width = 6, height = 6)
        baf.model.view(cellularity = cellularity, ploidy = ploidy,
            segs = seg.res[!segs.is.xy, ])
    dev.off()
    if (!is.null(cp.table)){
        alt.sol <- alternative.cp.solutions(cp.table)
        write.table(alt.sol, file = alt.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
        pdf(afit.file)
        for (sol in 1:nrow(alt.sol)){
            baf.model.view(cellularity = alt.sol$cellularity[sol],
                ploidy = alt.sol$ploidy[sol], segs = seg.res[!segs.is.xy, ])
        }
        dev.off()
    }
    file_conn <- file(log.file)
    writeLines(c(date(), paste("Sequenza version:",
        packageVersion("sequenza"), sep = " ")), file_conn)
    close(file_conn)
}
