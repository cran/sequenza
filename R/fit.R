sequenza.fit <- function(sequenza.extract, female = TRUE,
    N.ratio.filter = 10, N.BAF.filter = 1, segment.filter = 3e6,
    mufreq.treshold = 0.10, XY = c(X = "X", Y = "Y"),
    cellularity = seq(0.1, 1, 0.01), ploidy = seq(1, 7, 0.1),
    ratio.priority = FALSE, method = "baf",
    priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24,
    mc.cores = getOption("mc.cores", 2L)){

    if (is.null(chromosome.list)) {
        mut.all <- do.call(rbind, sequenza.extract$mutations)
        mut.all <- na.exclude(mut.all)
        segs.all <- do.call(rbind, sequenza.extract$segments)
    } else {
        mut.all <- do.call(rbind,
            sequenza.extract$mutations[chromosome.list])
        mut.all <- na.exclude(mut.all)
        segs.all <- do.call(rbind,
            sequenza.extract$segments[chromosome.list])
    }
    segs.len <- segs.all$end.pos - segs.all$start.pos
    avg.depth.ratio <- sequenza.extract$avg.depth.ratio

    if (method == "baf") {

        avg.sd.ratio <- sum(segs.all$sd.ratio * segs.all$N.ratio,
            na.rm = TRUE) / sum(segs.all$N.ratio, na.rm = TRUE)
        avg.sd.Bf <- sum(segs.all$sd.BAF * segs.all$N.BAF, na.rm = TRUE) /
            sum(segs.all$N.BAF, na.rm = TRUE)
        segs.all$sd.BAF[segs.all$sd.BAF == 0] <- max(
            segs.all$sd.BAF, na.rm = TRUE)
        segs.all$sd.ratio[segs.all$sd.ratio == 0] <- max(
            segs.all$sd.ratio, na.rm = TRUE)

        segs.filt <- segs.all$N.ratio > N.ratio.filter &
            segs.all$N.BAF > N.BAF.filter
        segs.filt <- segs.len >= segment.filter & segs.filt
        if (female){
            segs.is.xy <- segs.all$chromosome == XY["Y"]
        } else{
            segs.is.xy <- segs.all$chromosome %in% XY
        }
        filt.test  <- segs.filt & !segs.is.xy
        seg.test   <- segs.all[filt.test, ]
        seg.len.mb <- segs.len[filt.test] / 1e6
        baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio,
            sd.ratio = seg.test$sd.ratio, weight.ratio = seg.len.mb,
            sd.Bf = seg.test$sd.BAF, weight.Bf = seg.len.mb,
            avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
            ploidy = ploidy, priors.table = priors.table,
            mc.cores = mc.cores, ratio.priority = ratio.priority)
    } else if (method == "mufreq") {
        mut.filt <- mut.all$F >= mufreq.treshold
        if (female){
            mut.is.xy  <- mut.all$chromosome == XY["Y"]
        } else{
            mut.is.xy  <- mut.all$chromosome %in% XY
        }
        filt.test  <- mut.filt & !mut.is.xy
        mut.test   <- mut.all[filt.test, ]
        w.mufreq   <- round(mut.test$good.reads, 0)
        mufreq.model.fit(mufreq = mut.test$F,
            depth.ratio = mut.test$adjusted.ratio,
            weight.ratio = 2 * w.mufreq, weight.mufreq = w.mufreq,
            avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
            ploidy = ploidy, priors.table = priors.table, mc.cores = mc.cores)
    } else {
      stop("The only available methods are \"baf\" and \"mufreq\"")
   }
}
