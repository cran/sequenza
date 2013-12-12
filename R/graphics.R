cp.plot <- function (cp.table, ...) {
   #colorgram(x = cp.table$x, y = cp.table$y, z = log(cp.table$z),
   #    colFn = heat, map = map, outlier = outlier, las = 1,
   #    xlab = "ploidy", ylab = "Cellularity", zlab = "log-likelihood",
   #    ...)
   colorgram(x = cp.table$x, y = cp.table$y, z = matrix(rank(cp.table$z),
                                                        nrow = nrow(cp.table$z)), colFn = colorRampPalette(c('white', 'lightblue')),
             las = 1, xlab = "Ploidy", ylab = "Cellularity", zlab = "Rank likelihood",
             ...)
}

cp.plot.contours <- function(cp.table, likThresh = c(0.95),
                             col = palette(), legend.pos = 'bottomright', pch = 18, ...) {
   znormsort <- sort(cp.table$z, decreasing = TRUE)
   znormcumLik <- cumsum(znormsort)
   n <- sapply(likThresh, function(x) sum(znormcumLik < x) + 1)
   LikThresh <- znormsort[n]
   names(LikThresh) <- paste0(likThresh * 100, '%')

   contour(cp.table, levels = znormsort[n], col = col,
           drawlabels = FALSE,
           xlab= "Ploidy", ylab = "Cellularity", ...)
   max.xy <- which(cp.table$z == max(cp.table$z), arr.ind = TRUE)
   points(x = cp.table$x[max.xy[, "row"]],
          y = cp.table$y[max.xy[, "col"]], pch = pch)
   if(!is.na(legend.pos)) {
      legend(legend.pos, legend = c(paste("C.R.", names(LikThresh), sep = " "), "Point estimate"),
             col = c(col[1:length(LikThresh)], "black"), lty = c(rep(1, length(LikThresh)), NA),
             pch = c(rep(NA, length(LikThresh)), pch), border = NA, bty = "n")
   }
   invisible(LikThresh)
}

# plot.fit.model <- function(mufreq.tab, cellularity, ploidy, chr23 = "XY",
#                            cn.ratio.range = c(0.5:2), avg.depth.ratio = avg.depth.ratio,
#                            cex.m = 1, cex.d = 1, ...) {
#    xy.index   <- mufreq.tab$chr == "chrX" | mufreq.tab$chr == "chrY"
#    plot(x = mufreq.tab$F, y = mufreq.tab$adjusted.ratio,
#         xlab = "mutation frequency", ylab = "depth.ratio",
#         las = 1, type="n", ...)
#    points(x = mufreq.tab$F[!xy.index], y = mufreq.tab$adjusted.ratio[!xy.index],
#           pch = 19, col = "blue", cex = cex.d)
#    types      <- types.matrix(cn.ratio.range = cn.ratio.range, chr23 = chr23)
#    if (length(which(xy.index)) >= 1) {
#       if (chr23 == "XY") {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "green")
#       } else {
#          points(x = mufreq.tab$F[xy.index], y = mufreq.tab$adjusted.ratio[xy.index], pch = 19, col = "blue")
#       }
#    }
#    points.fit <-model.points(cellularity = cellularity, ploidy = ploidy,
#                              types = types, avg.depth.t = avg.depth.t,
#                              avg.depth.r = avg.depth.r)
#    points(points.fit,pch = 19, col = "red", cex = cex.m)
#    if (length(which(xy.index)) >= 1 & chr23 == "XY") {
#       legend("bottomright", c("Male chr X/Y",paste(paste("C", cellularity,sep = " : "), paste("P", ploidy,sep = " : "), sep = "; ")), pch = 19, col = c("green","red"))
#    } else {
#       legend("bottomright", paste(paste("C", cellularity,sep = " : "), paste("P", ploidy,sep = " : "), sep = "; "), pch = 19, col = "red")
#    }
# }

# plot.gc.depth <- function(d.values, gc.contents, ...) {
#    gc.values <- sort(unique(gc.contents))
#    d.list <- list()
#    #sizevect <- rep(0,length(gc.values))
#    for ( i in 1:length(gc.values)) {
#       d.list[[i]] <- d.values[gc.contents == gc.values[i]]
#       #sizevect[i]    <- length(d.list[[i]])
#    }
#    bxplot(d.list, names = gc.values, ...)
# }

plotWindows <- function(abf.window, m.lty = 1, m.lwd = 3,
                         m.col = "black", q.bg = "lightblue", log2.plot = FALSE,
                         n.min = 1, xlim, ylim, add = FALSE, ...) {
   if (log2.plot == TRUE) {
      abf.window[, c(3, 4, 5)] <- log2(abf.window[, c(3, 4, 5)])
   }
   if(!add) {
      if(missing(xlim))
         xlim <- c(abf.window[1, 1], abf.window[nrow(abf.window), 2])
      if(missing(ylim))
         ylim <- c(min(abf.window[, 4], na.rm = TRUE), max(abf.window[, 5], na.rm = TRUE))
      plot(xlim, ylim, type = "n", ...)
   }
   abf.window <- abf.window[abf.window[, 6] >= n.min, ]
   rect(xleft = abf.window[, 1], ybottom = abf.window[, 4],
        xright = abf.window[, 2], ytop = abf.window[, 5],
        col = q.bg, border = NA)
   segments(y0 = abf.window[, 3], x0 = abf.window[, 1] , x1 = abf.window[, 2], lty = m.lty, lwd = m.lwd, col = m.col)

}

chromosome.view <- function(baf.windows, ratio.windows, mut.tab = NULL, segments = NULL,  min.N.baf = 1, min.N.ratio = 1e4,
                            main = "", vlines = FALSE, legend.inset = c(-20 * strwidth("a", units = 'figure'), 0), BAF.style = "none",
                            CNn = 2, cellularity = NULL, ploidy = NULL, avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24", model.col = 1,
                            x.chr.space = 10) {
   make.polygons <- function(segments, model.baf) {
      max.B      <- max(model.baf$B[model.baf$CNt == max(segments$CNt)])
      mat.polygs <- matrix(ncol = max.B+1, nrow = nrow(segments))
      colnames(mat.polygs) <- 0:max.B
      get.B <- function (CNt, B) model.baf$BAF[model.baf$CNt == CNt & model.baf$B == B]
      polyg.coords <- sapply( 0:max.B, FUN = function (k) as.numeric(sapply(segments$CNt, FUN = function(i) get.B(i, k))))
      polyg.coords[is.na(polyg.coords)] <- 1
      polyg.coords <- cbind(polyg.coords, 1)
      polyg.pos    <- segments[, c("start.pos", "end.pos")]
      edge1        <- polyg.pos$end.pos[-nrow(polyg.pos)]
      edge2        <- polyg.pos$start.pos[-1]
      no.dat       <- c(1:nrow(polyg.pos))[edge2 - edge1 >= 1e6]
      v.gaps       <- apply(rbind(edge1[no.dat], edge2[no.dat]), 2, mean)
      v.gaps       <- cbind(start.pos = v.gaps, end.pos = v.gaps)
      polyg.p.new  <- rbind(polyg.pos, v.gaps)
      polyg.c.new  <- polyg.coords
      for (i in no.dat) {
         polyg.c.new  <- rbind(polyg.c.new, 0)
      }
      new_idx      <- c(seq_along(polyg.pos$start.pos), no.dat+0.5)
      polyg.pos    <- polyg.p.new[order(new_idx), ]
      polyg.coords <- polyg.c.new[order(new_idx), ]
      Xs      <- unlist(lapply(1:nrow(polyg.coords), function(k) polyg.pos[k, ]))
      #color <- colorRampPalette(c("grey99", "grey20"))( max.B + 1 )
      color <- gray.colors((max.B + 1), start = 0.5, end = 0.9, alpha = 0.3)
      extra.x <- c(max(Xs), min(Xs))
      extra.y = c(0 ,0)
      for (k in 1:ncol(polyg.coords)) {
         Ys <- unlist(lapply(polyg.coords[, k], function (i) rep(i, 2)))
         polygon(x = c(Xs,extra.x), y = c(Ys,extra.y), col = color[k], border = NA)
         extra.y <- rev(Ys)
         extra.x <- rev(Xs)
      }
   }
   if (is.null(segments)) {
      data.model <- NULL
   } else {
      if ("CNt" %in% colnames(segments)) {
         if (length(c(cellularity, ploidy, avg.depth.ratio)) != 3) {
            data.model <- NULL
         } else {
            data.model     <- list()
            CNt.max        <- max(segments$CNt, na.rm = TRUE) + 1
            CNt.min        <- 0
            data.model$baf <- theoretical.baf(CNn = CNn, CNt = CNt.max, cellularity = cellularity)
            if (CNn == 2) {
               data.model$baf <- rbind(c(0,0,0.5,0), data.model$baf)
            } else {
               data.model$baf <- rbind(c(0,0,1,0), data.model$baf)
            }
            types          <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNn = CNn)
            data.model$muf <- cbind(types, model.points(cellularity = cellularity, ploidy = ploidy,
                                                   types = types, avg.depth.ratio = avg.depth.ratio))
         }
      } else {
         data.model <- NULL
      }
   }
   if (is.null(mut.tab)) {
      par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0), mfcol = c(2,1), xaxt='n')
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
      xlim <- c(min.x, max.x)
   } else {
      min.x <- min(c(min(baf.windows$start), min(ratio.windows$start), min(mut.tab$n.base)))
      max.x <- max(c(max(baf.windows$end), max(ratio.windows$end), max(mut.tab$n.base)))
      xlim <- c(min.x, max.x)
      par(mar = c(0, 4, 0, 10), oma = c(5, 0, 4, 0), mfcol = c(3,1), xaxt='n', xpd = TRUE)
      mutation.colors <- c(
         'A>C' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'T>G' = rgb(red =   0, green = 178, blue = 238, alpha = 120, maxColorValue = 255),
         'A>G' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'T>C' = rgb(red = 255, green =  64, blue =  64, alpha = 120, maxColorValue = 255),
         'A>T' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'T>A' = rgb(red =  34, green = 139, blue =  34, alpha = 120, maxColorValue = 255),
         'C>A' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'G>T' = rgb(red = 139, green =  90, blue =   0, alpha = 120, maxColorValue = 255),
         'C>G' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'G>C' = rgb(red = 127, green =   0, blue = 255, alpha = 120, maxColorValue = 255),
         'C>T' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255),
         'G>A' = rgb(red = 255, green = 215, blue =   0, alpha = 120, maxColorValue = 255)
      )
      plot(x = mut.tab$n.base, y = mut.tab$F,
           ylab = "Mutant allele frequency", las = 1, pch = 19,
           col = c(mutation.colors, 'NA' = NA)[as.character(mut.tab$mutation)],
           ylim = c(min(mut.tab$F, na.rm = TRUE), 1), xlim = xlim)
      unique.colors <- unique(mutation.colors)
      labels <- sapply(unique.colors, function(a) paste(names(mutation.colors)[mutation.colors == a], collapse = ", "))
      #legend("topleft", legend = labels, fill = unique.colors, border = NA, bty = "n")
      legend(y = "center", x  = "right", legend = labels,
             inset = legend.inset,
             fill = unique.colors, border = NA, bty = "n")
      if (!is.null(segments)){
         if (vlines) {
            abline(v = segments$end.pos, lwd = 1, lty = 2)
         }
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
                segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                         y0 = unique(data.model$muf$mufreqs[data.model$muf$CNt == segments$CNt[i]]), lwd = model.lwd, lty = model.lty, col = model.col)
            }
         }
      }

   }
   if (!is.null(segments)){
      plot(ylab = "B allele frequency", type = "n",
           x = xlim, y = c(0, 0.5), las = 1)
      if (BAF.style == "blocks") {
         if (!is.null(data.model)) {
            make.polygons(segments, data.model$baf)
            axis(side = 4, line = 0, las = 1,
                 labels = data.model$baf$B[data.model$baf$CNt == segments$CNt[nrow(segments)]],
                 at = data.model$baf$BAF[data.model$baf$CNt == segments$CNt[nrow(segments)]])
            mtext(text = "Number of B alleles", side = 4, line = 2, cex = par("cex.lab")*par("cex"))
         }
      }
      plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf, add = TRUE)
      if (vlines) {
         abline(v = segments$end.pos, lwd = 1, lty = 2)
      }
      segments(x0 = segments$start.pos, y0 = segments$Bf, x1=segments$end.pos, y1 = segments$Bf, col = "red", lwd = 3)
      if (BAF.style == "lines") {
         if (!is.null(data.model)) {
            for (i in 1:nrow(segments)) {
               segments(x0 = segments$start.pos[i], x1 = segments$end.pos[i],
                        y0 = unique(data.model$baf$BAF[data.model$baf$CNt == segments$CNt[i]]), lwd = model.lwd, lty = model.lty, col = model.col)
            }
         }
      }
   }
   else {
      plotWindows(baf.windows, ylab = "B allele frequency",
                  xlim = xlim, ylim = c(0, 0.5), las = 1,
                  n.min = min.N.baf)
   }
   plotWindows(ratio.windows, ylab = "Depth ratio",
               las = 1, n.min = min.N.ratio, ylim = c(0, 2.5))
   if (!is.null(segments)){
      if (vlines) {
         abline(v = segments$end.pos, lwd = 1, lty = 2)
      }
      segments(x0 = segments$start.pos, y0 = segments$depth.ratio, x1=segments$end.pos, y1 = segments$depth.ratio, col = "red", lwd = 3)
      if (!is.null(data.model)) {
         ratios.theoric <- unique(data.model$muf[,c('CNt', 'depth.ratio')])

         segments(x0 = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
                  x1 = rep(max(segments$end.pos, na.rm = TRUE), times = nrow(ratios.theoric)),
                  y0 = ratios.theoric$depth.ratio, lwd = model.lwd, lty = model.lty, col = model.col)
         #text(x = rep(min(segments$start.pos, na.rm =TRUE), times = nrow(ratios.theoric)),
         #     y = ratios.theoric$depth.ratio, labels = ratios.theoric$CNt, pos = 2, offset = 0.5, cex = 0.8)
         axis(labels = as.character(ratios.theoric$CNt), side = 4, line = 0, las = 1,
              at = ratios.theoric$depth.ratio)
         mtext(text = "Copy number", side = 4, line = 2, cex = par("cex.lab")*par("cex"))
      }
   }
   par(xaxt='s')
   axis(labels = as.character(round(seq(xlim[1]/1e6, xlim[2]/1e6, by = x.chr.space), 0)), side = 1 , line = 0,
        at = seq(xlim[1], xlim[2], by = 1e6 * x.chr.space), outer = FALSE, cex = par("cex.axis")*par("cex"))
   mtext("Position (Mb)", side = 1, line = 3, outer = FALSE, cex = par("cex.lab")*par("cex"))
   mtext(main, 3, outer = TRUE, cex = par("cex.main")*par("cex"), line = 2)
}

#genome.view <- function(baf.windows, ratio.windows, segments = NULL, main = "",
#                            min.N.baf = 1, min.N.ratio = 1e4, CNn = rep(2, length(ratio.windows)),
#                            cellularity = NULL, ploidy = NULL, avg.depth.ratio = NULL) {
#   chr.metrics <- list()
#   for (i in 1:length(ratio.windows)) {
#      chr.metrics[[i]] <- range(ratio.windows[[i]]$mean, na.rm = TRUE)
#   }
#   chr.metrics <- do.call(rbind, chr.metrics)
#   x0 <- chr.metrics[1,1]
#
#}

genome.view <- function(seg.cn, info.type = "AB", ...) {
   chr.order <- unique(seg.cn$chromosome)
   seg.list  <- split(x = seg.cn[,c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")],
                      f = seg.cn$chromosome)
   seg.list  <- seg.list[order(order(chr.order))]
   seg.max   <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
   seg.pos   <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
   seg.max   <- cumsum(as.numeric(do.call(rbind, seg.max)))
   chr.offset <- 0
   for (i in 1:length(seg.pos)){
      seg.pos[[i]] <- seg.pos[[i]] + chr.offset
      colnames(seg.pos[[i]]) <- c("abs.start","abs.end")
      chr.offset   <- seg.max[i]
   }
   seg.max      <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
   abs.list     <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
   abs.segments <- do.call(rbind, abs.list)
   if (info.type == "AB") {
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(-0.1, (max(abs.segments$A)+0.1)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n',  yaxt='n', xaxs = "i", ...)
      axis(labels = 0:max(abs.segments$A),
           at = 0:max(abs.segments$A),
           side = 2, line = 0, las = 1)
      #abline(h = c(0:max(abs.segments$A)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$B-0.1), y1 = (abs.segments$B-0.1), col="blue", lwd = 5, lend = 1)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = (abs.segments$A+0.1), y1 = (abs.segments$A+0.1), col="red", lwd = 5, lend = 1)
   } else {
      plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
           y = c(min(abs.segments$CNt), max(abs.segments$CNt)), type = "n",
           ylab = "Copy number", xlab = "Position (Mb)",
           xaxt='n', yaxt = 'n', xaxs = "i", ...)
      axis(labels = min(abs.segments$CNt):max(abs.segments$CNt),
           at = min(abs.segments$CNt):max(abs.segments$CNt),
           side = 2, line = 0, las = 1)
      #abline(h = c(min(abs.segments$CNt):max(abs.segments$CNt)), lty = 2)
      segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
               y0 = abs.segments$CNt, y1= abs.segments$CNt, col="red", lwd = 5, lend = 1)
   }
   abline(v = c(0, seg.max), lty = 3)
   for (i in 1:length(abs.list)){
      max.pos <- nrow(abs.list[[i]])
      mtext(chr.order[i], side = 3, line = 0,
            at = sum(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos])/2)
      #axis(labels = as.character(round(seq(abs.list[[i]]$start.pos[1]/1e6, abs.list[[i]]$end.pos[max.pos]/1e6, by = 20), 0)),
      #     at = seq(abs.list[[i]]$abs.start[1], abs.list[[i]]$abs.end[max.pos], by = 2e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
      #     side = 1 , line = 0)
   }
   axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1]/1e6, abs.list[[1]]$end.pos[nrow(abs.list[[1]])]/1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1], abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7), outer = FALSE, cex = par("cex.axis")*par("cex"),
        side = 1 , line = 1)
}
