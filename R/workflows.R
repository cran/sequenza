sequenza.extract <- function(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80, kmin = 10,
                             mufreq.treshold = 0.10, min.reads = 40, max.mut.types = 1,
                             min.type.freq = 0.9){
   gc.stats <- gc.sample.stats(file, gz = gz)
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
      if (nrow(abf.het) > 0) {
         abf.b.win <- windowValues(x = abf.het$Bf,
                                   positions = abf.het$n.base,
                                   chromosomes = abf.het$chromosome,
                                   window = window, overlap = overlap,
                                   weight = abf.het$good.s.reads)
         breaks    <- NULL
         breaks    <- try(find.breaks(abf.het, gamma = gamma, 
                                      kmin = kmin, baf.thres = c(0, 0.5)),
                          silent = FALSE)
         if (!is.null(breaks)){
            seg.s1    <- segment.breaks(abf.data, breaks = breaks)            
         } else {
            seg.s1 <- segment.breaks(abf.data,
                                     breaks = data.frame(chrom = chr,
                                                         start.pos = min(abf.data$n.base, na.rm = TRUE),
                                                         end.pos = max(abf.data$n.base, na.rm = TRUE)))
         }
         
      } else {
         abf.b.win <- list()
         abf.b.win[[1]] <- data.frame(start = min(abf.data$n.base, na.rm = TRUE),
                                      end = max(abf.data$n.base, na.rm = TRUE), mean = 0.5,
                                      q0 = 0.5,  q1 = 0.5, N = 1)
         seg.s1 <- segment.breaks(abf.data,
                                  breaks = data.frame(chrom = chr,
                                                      start.pos = min(abf.data$n.base, na.rm = TRUE),
                                                      end.pos = max(abf.data$n.base, na.rm = TRUE)))
                                  
      }
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
                         cellularity = seq(0.1,1,0.01), ploidy = seq(1, 7, 0.1), ratio.priority = FALSE,
                         priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24,
                         mc.cores = getOption("mc.cores", 2L)){
   if (is.null(chromosome.list)) {
      segs.all      = do.call(rbind, sequenza.extract$segments)
   } else {
      segs.all      = do.call(rbind, sequenza.extract$segments[chromosome.list])      
   }
   #mut.all       = do.call(rbind, sequenza.extract$mutations)
   #mut.all       = na.exclude(mut.all)
   segs.len      = segs.all$end.pos - segs.all$start.pos
   segs.filt     = segs.len >= segment.filter
   avg.depth.ratio = mean(sequenza.extract$gc$adj[,2])
   if (female == FALSE){
      segs.is.xy = segs.all$chromosome == XY["X"] | segs.all$chromosome == XY["Y"]
      #mut.is.xy  = mut.all$chromosome == XY["X"] | mut.all$chromosome == XY["Y"]
   } else{
      segs.is.xy = segs.all$chromosome == XY["Y"]
      #mut.is.xy  = mut.all$chromosome == XY["Y"]
   }
   filt.test  = segs.filt & !segs.is.xy
   seg.test   = segs.all[filt.test, ]
   weights.seg = round(segs.len[filt.test] / 1e6, 0) + 150
   baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio,
                 weight.ratio = 2 * weights.seg, weight.Bf = weights.seg,
                 avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
                 ploidy = ploidy, priors.table = priors.table,
                 mc.cores = mc.cores, ratio.priority = ratio.priority)
}

sequenza.results <- function(sequenza.extract, sequenza.fit = NULL, sample.id, out.dir = './',
                             cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                             ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
   cp.file   <- paste(sample.id, "CP_contours.pdf", sep = '_')
   cint.file <- paste(sample.id, "confints_CP.txt", sep = '_')
   chrw.file <- paste(sample.id, "chromosome_view.pdf", sep = '_')
   geno.file <- paste(sample.id, "genome_view.pdf", sep = '_')
   cn.file   <- paste(sample.id, "CN_bars.pdf", sep = '_')
   muts.file <- paste(sample.id, "mutations.txt", sep = '_')
   segs.file <- paste(sample.id, "segments.txt", sep = '_')
   robj.extr <- paste(sample.id, "sequenza_extract.RData", sep= '_')
   robj.fit  <- paste(sample.id, "sequenza_fit.RData", sep= '_')  
   avg.depth.ratio <- mean(sequenza.extract$gc$adj[, 2])
   dir.create(path = out.dir, showWarnings = FALSE,
              recursive = TRUE)
   assign(x = paste0(sample.id,"_sequenza_extract"), value = sequenza.extract)
   save(list = paste0(sample.id,"_sequenza_extract"), file = paste(out.dir, robj.extr, sep = "/")) 
   if (is.null(sequenza.fit) & (is.null(cellularity) | is.null(ploidy))){
      stop("Either the sequenza.fit or both cellularity and ploidy argument are required.")
   }
   if (!is.null(sequenza.fit)){
      assign(x = paste0(sample.id,"_sequenza_fit"), value = sequenza.fit)
      save(list = paste0(sample.id,"_sequenza_fit"), file = paste(out.dir, robj.fit, sep = "/"))       
      cint <- get.ci(sequenza.fit)
      pdf(paste(out.dir, cp.file, sep = "/"))
         cp.plot(sequenza.fit)
         cp.plot.contours(sequenza.fit, add=T, likThresh = c(0.95), col="red", pch = 20)
         if (!is.null(cellularity) | !is.null(ploidy)) {
            if (is.null(cellularity)) cellularity <- cint$max.y
            if (is.null(ploidy)) ploidy <- cint$max.x
            points(x = ploidy, y = cellularity , pch = 5)
            text(x = ploidy, y = cellularity, 
                 labels = "User selection",
                 pos = 3, offset = 0.5)
         } else {
            cellularity <- cint$max.y
            ploidy <- cint$max.x
         }
      dev.off()
   }
   seg.tab     <- na.exclude(do.call(rbind, sequenza.extract$segments[chromosome.list = 1:24]))
   mut.tab     <- na.exclude(do.call(rbind, sequenza.extract$mutations[chromosome.list = 1:24]))
   if (female == FALSE){
      segs.is.xy = seg.tab$chromosome == XY["X"] | seg.tab$chromosome == XY["Y"]
      mut.is.xy  = mut.tab$chromosome == XY["X"] | mut.tab$chromosome == XY["Y"]
   } else{
      segs.is.xy = seg.tab$chromosome == XY["Y"]
      mut.is.xy  = mut.tab$chromosome == XY["Y"]
   }

   cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
                            depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio,
                            ratio.priority = ratio.priority, CNn = 2)
   seg.res     <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
   if (female == FALSE){
      cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[segs.is.xy], CNt.max = CNt.max,
                               depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                               cellularity = cellularity, ploidy = ploidy,
                               avg.depth.ratio = avg.depth.ratio,
                               ratio.priority = TRUE, CNn = 1)
      seg.xy     <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
      seg.res    <- rbind(seg.res, seg.xy)     
   }
   write.table(seg.res, paste(out.dir, segs.file, sep = "/"),
               col.names = TRUE, row.names = FALSE, sep="\t")   

   mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy], CNt.max = CNt.max,
                            depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio, CNn = 2)
   mut.res     <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
   if (female == FALSE){
      mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy], CNt.max = CNt.max,
                                  depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                                  cellularity = cellularity, ploidy = ploidy,
                                  avg.depth.ratio = avg.depth.ratio, CNn = 1)
      mut.xy     <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
      mut.res    <- rbind(mut.res, mut.xy)     
   }
   write.table(mut.res, paste(out.dir, muts.file, sep = "/"),
               col.names = TRUE, row.names = FALSE, sep="\t")   
   
   pdf(paste(out.dir, chrw.file, sep = "/"))
   for (i in unique(seg.res$chromosome)) {
      
      if (female == FALSE & (i == XY["X"] | i == XY["Y"])){
         CNn = 1
      } else {
         CNn <- 2
      }
      chromosome.view(mut.tab = sequenza.extract$mutations[[i]],
                      baf.windows = sequenza.extract$BAF[[i]],
                      ratio.windows = sequenza.extract$ratio[[i]], BAF.style="lines",
                      cellularity = cellularity, ploidy = ploidy, main = i,
                      segments = seg.res[seg.res$chromosome == i, ],
                      avg.depth.ratio = avg.depth.ratio, CNn = CNn, min.N.ratio = 1)
   }
   dev.off()
   pdf(paste(out.dir, geno.file, sep = "/"))
      genome.view(seg.res)
      genome.view(seg.res, "CN")
   dev.off()
   barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos,
                     CNt = seg.res$CNt)
   cn.sizes <- split(barscn$size,barscn$CNt)
   cn.sizes <- sapply(cn.sizes, 'sum')
   pdf(paste(out.dir, cn.file, sep = "/"))
   barplot(round(cn.sizes/sum(cn.sizes)* 100, 0), names = names(cn.sizes), las = 1,
                      ylab = "Percentage (%)", xlab = "copy number")
   dev.off()
   
   ## Write down the results.... ploidy etc...
   if (!is.null(sequenza.fit)){
      res.tab = data.frame(cellularity     = c(cint$confint.y[1], cint$max.y[1], cint$confint.y[2]),
                           ploidy.estimate = c(cint$confint.x[1], cint$max.x[1], cint$confint.x[2]),
                           ploidy.mean.cn  = weighted.mean(x = as.integer(names(cn.sizes)), w = cn.sizes))
      write.table(res.tab, paste(out.dir, cint.file, sep = "/"), col.names = TRUE,
                  row.names = FALSE, sep = "\t")
   }
}
