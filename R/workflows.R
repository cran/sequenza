sequenza.extract <- function(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80, kmin = 10,
                             mufreq.treshold = 0.10, min.reads = 40, max.mut.types = 1,
                             min.type.freq = 0.9){
   gc.stats <- gc.sample.stats(file)
   chr.vect <- as.character(gc.stats$file.metrics$chr)
   gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
   windows.baf   <- list()
   windows.ratio <- list()
   mutation.list <- list()
   segments.list <- list()
   for (chr in chr.vect){
      file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
      abf.data   <- read.abfreq(file, gz = gz, n.lines = c(file.lines$start, file.lines$end))
      abf.data$adjusted.ratio <- round(abf.data$depth.ratio / gc.vect[as.character(abf.data$GC.percent)], 3)
      abf.hom <- abf.data$ref.zygosity == 'hom'
      abf.het <- abf.data[!abf.hom, ]
      abf.r.win <- windowValues(x = abf.data$adjusted.ratio,
                                positions = abf.data$n.base,
                                chromosomes = abf.data$chromosome,
                                window = window, overlap = overlap,
                                weight = abf.data$depth.normal)
      abf.b.win <- windowValues(x = abf.het$Bf,
                                positions = abf.het$n.base,
                                chromosomes = abf.het$chromosome,
                                window = window, overlap = overlap,
                                weight = abf.het$good.s.reads)
      breaks    <- find.breaks(abf.het, gamma = gamma, kmin = kmin, baf.thres = c(0, 0.5))
      seg.s1    <- segment.breaks(abf.data, breaks = breaks)
      mut.tab   <- mutation.table(abf.data, mufreq.treshold = mufreq.treshold,
                                  min.reads = min.reads, max.mut.types = max.mut.types,
                                  min.type.freq = min.type.freq, segments = seg.s1)
      windows.baf[[which(chr.vect == chr)]]   = abf.b.win[[1]]
      windows.ratio[[which(chr.vect == chr)]] = abf.r.win[[1]]
      mutation.list[[which(chr.vect == chr)]] = mut.tab
      segments.list[[which(chr.vect == chr)]] = seg.s1
   }
   names(windows.baf)   <- chr.vect
   names(windows.ratio) <- chr.vect
   names(mutation.list) <- chr.vect
   names(segments.list) <- chr.vect
   return(list(BAF = windows.baf, ratio = windows.ratio, mutations = mutation.list,
               segments = segments.list, chromosomes = chr.vect, gc = gc.stats))
}

sequenza.fit <- function(sequenza.extract, female = TRUE, segment.filter = 1e7, XY = c(X = "X", Y = "Y"),
                         cellularity = seq(0.1,1,0.01), ploidy = seq(1, 7, 0.1),
                         priors.table = data.frame(CN = 2, value = 2),
                         mc.cores = getOption("mc.cores", 2L)){
   segs.all      = do.call(rbind, sequenza.extract$segments)
   mut.all       = do.call(rbind, sequenza.extract$mutations)
   mut.all       = na.exclude(mut.all)
   segs.len      = segs.all$end.pos - segs.all$start.pos
   segs.filt     = segs.len >= segment.filter
   avg.depth.ratio = mean(sequenza.extract$gc$adj[,2])
   if (female == FALSE){
      segs.is.xy = segs.all$chromosome == XY["X"] | segs.all$chromosome == XY["Y"]
      mut.is.xy  = mut.all$chromosome == XY["X"] | mut.all$chromosome == XY["Y"]
   } else{
      segs.is.xy = segs.all$chromosome == XY["Y"]
      mut.is.xy  = mut.all$chromosome == XY["Y"]
   }
   filt.test  = segs.filt & !segs.is.xy
   seg.test   = segs.all[filt.test, ]
   weights.seg = round(segs.len[filt.test] / 1e6, 0) + 150
   baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio,
                 weight.ratio = 2 * weights.seg, weight.Bf = weights.seg,
                 avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
                 ploidy = ploidy, priors.table = priors.table,
                 mc.cores = mc.cores)
}
