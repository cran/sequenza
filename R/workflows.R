sequenza.extract <- function(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80, kmin = 10,
                             mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
                             min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
                             min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
                             weighted.mean = TRUE){
   gc.stats <- gc.sample.stats(file, gz = gz)
   chr.vect <- as.character(gc.stats$file.metrics$chr)
   gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
   windows.baf   <- list()
   windows.ratio <- list()
   mutation.list <- list()
   segments.list <- list()
   coverage.list <- list()
   if (is.null(chromosome.list)) {
     chromosome.list <- chr.vect
   } else {
     chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
   }
   for (chr in chromosome.list){
      if (verbose){
         message("Processing ", chr, ": ", appendLF = FALSE) 
      }
      file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
      seqz.data   <- read.seqz(file, gz = gz, n.lines = c(file.lines$start, file.lines$end))
      seqz.data$adjusted.ratio <- round(seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)], 3)
      seqz.hom <- seqz.data$zygosity.normal == 'hom'
      seqz.het <- seqz.data[!seqz.hom, ]
      het.filt <- seqz.het$good.reads >= min.reads.baf
      seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
                                positions = seqz.data$position,
                                chromosomes = seqz.data$chromosome,
                                window = window, overlap = overlap,
                                weight = seqz.data$depth.normal)
      if (nrow(seqz.het) > 0) {
         seqz.b.win <- windowValues(x = seqz.het$Bf[het.filt],
                                   positions = seqz.het$position[het.filt],
                                   chromosomes = seqz.het$chromosome[het.filt],
                                   window = window, overlap = overlap,
                                   weight = seqz.het$good.reads[het.filt])
         breaks    <- NULL
         breaks    <- try(find.breaks(seqz.het, gamma = gamma, 
                                      kmin = kmin, baf.thres = c(0, 0.5)),
                          silent = FALSE)
         if (!is.null(breaks)){
            seg.s1    <- segment.breaks(seqz.tab = seqz.data, breaks = breaks,
                                        min.reads.baf = min.reads.baf, weighted.mean = weighted.mean)            
         } else {
            seg.s1 <- segment.breaks(seqz.data,
                                     breaks = data.frame(chrom = chr,
                                                         start.pos = min(seqz.data$position, na.rm = TRUE),
                                                         end.pos = max(seqz.data$position, na.rm = TRUE)),
                                     weighted.mean = weighted.mean)
         }
         
      } else {
         seqz.b.win <- list()
         seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position, na.rm = TRUE),
                                      end = max(seqz.data$position, na.rm = TRUE), mean = 0.5,
                                      q0 = 0.5,  q1 = 0.5, N = 1)
         seg.s1 <- segment.breaks(seqz.data,
                                  breaks = data.frame(chrom = chr,
                                                      start.pos = min(seqz.data$position, na.rm = TRUE),
                                                      end.pos = max(seqz.data$position, na.rm = TRUE)),
                                  weighted.mean = weighted.mean)
                                  
      }
      mut.tab   <- mutation.table(seqz.data, mufreq.treshold = mufreq.treshold,
                                  min.reads = min.reads, min.reads.normal = min.reads.normal,
                                  max.mut.types = max.mut.types, min.type.freq = min.type.freq,
                                  min.fw.freq = min.fw.freq, segments = seg.s1)
      windows.baf[[which(chromosome.list == chr)]]   <- seqz.b.win[[1]]
      windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
      mutation.list[[which(chromosome.list == chr)]] <- mut.tab
      segments.list[[which(chromosome.list == chr)]] <- seg.s1
      coverage.list[[which(chromosome.list == chr)]] <- data.frame(sum = sum(as.numeric(seqz.data$depth.tumor),
                                                                       na.rm = TRUE),
                                                                 N  = length(seqz.data$depth.tumor) )
      if (verbose){
        
        message(nrow(mut.tab), ' variant calls; ',
                 nrow(seqz.het), ' heterozygous positions; ',
                 sum(seqz.hom), ' homozygous positions.') 
      }
   }
   names(windows.baf)   <- chromosome.list
   names(windows.ratio) <- chromosome.list
   names(mutation.list) <- chromosome.list
   names(segments.list) <- chromosome.list
   coverage.list <- do.call(rbind, coverage.list)
   coverage <- sum(coverage.list$sum) / sum(coverage.list$N)
   return(list(BAF = windows.baf, ratio = windows.ratio, mutations = mutation.list,
               segments = segments.list, chromosomes = chromosome.list, gc = gc.stats,
               avg.depth = round(coverage,0)))
}

sequenza.fit <- function(sequenza.extract, female = TRUE, segment.filter = 1e7, mufreq.treshold = 0.10, 
                         XY = c(X = "X", Y = "Y"), cellularity = seq(0.1, 1, 0.01), ploidy = seq(1, 7, 0.1),
                         ratio.priority = FALSE, method = "baf", priors.table = data.frame(CN = 2, value = 2),
                         chromosome.list = 1:24, mc.cores = getOption("mc.cores", 2L)){
   if (method == "baf") {
      if (is.null(chromosome.list)) {
         segs.all      <- do.call(rbind, sequenza.extract$segments)
      } else {
         segs.all      <- do.call(rbind, sequenza.extract$segments[chromosome.list])      
      }
      segs.len      <- segs.all$end.pos - segs.all$start.pos
      segs.filt     <- segs.len >= segment.filter
      avg.depth.ratio = mean(sequenza.extract$gc$adj[,2])
      if (female){
         segs.is.xy <- segs.all$chromosome == XY["Y"]
      } else{
         segs.is.xy <- segs.all$chromosome %in% XY
      }
      filt.test  <- segs.filt & !segs.is.xy
      seg.test   <- segs.all[filt.test, ]
      weights.seg <- round(segs.len[filt.test] / 1e6, 0) + 150
      baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio,
                    weight.ratio = 2 * weights.seg, weight.Bf = weights.seg,
                    avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
                    ploidy = ploidy, priors.table = priors.table,
                    mc.cores = mc.cores, ratio.priority = ratio.priority)
   } else if (method == "mufreq") {
      if (is.null(chromosome.list)) {
         mut.all       <- do.call(rbind, sequenza.extract$mutations)
         mut.all       <- na.exclude(mut.all)
      } else {
         mut.all       <- do.call(rbind, sequenza.extract$mutations[chromosome.list])
         mut.all       <- na.exclude(mut.all)
      }
      mut.filt     <- mut.all$F >= mufreq.treshold
      avg.depth.ratio = mean(sequenza.extract$gc$adj[,2])
      if (female){
         mut.is.xy  <- mut.all$chromosome == XY["Y"]
      } else{
         mut.is.xy  <- mut.all$chromosome %in% XY
      }
      filt.test  <- mut.filt & !mut.is.xy
      mut.test   <- mut.all[filt.test, ]
      w.mufreq   <- round(mut.test$good.reads, 0)
      mufreq.model.fit(mufreq = mut.test$F, depth.ratio = mut.test$adjusted.ratio,
                    weight.ratio = 2 * w.mufreq, weight.mufreq = w.mufreq,
                    avg.depth.ratio = avg.depth.ratio, cellularity = cellularity,
                    ploidy = ploidy, priors.table = priors.table,
                    mc.cores = mc.cores)    
   } else {
      stop("The only available methods are \"baf\" and \"mufreq\"")
   }
}

sequenza.results <- function(sequenza.extract, cp.table = NULL, sample.id, out.dir = getwd(),
                             cellularity = NULL, ploidy = NULL, female = TRUE, CNt.max = 20,
                             ratio.priority = FALSE, XY = c(X = "X", Y = "Y"), chromosome.list = 1:24){
   if(!file.exists(out.dir)) {
     dir.ok <- dir.create(path = out.dir, recursive = TRUE)
     if(!dir.ok) stop('Directory does not exist and cannot be created: ', out.dir)
   }
   makeFilename <- function(x) file.path(out.dir, paste(sample.id, x, sep = '_'))
   cp.file   <- makeFilename("CP_contours.pdf")
   cint.file <- makeFilename("confints_CP.txt")
   chrw.file <- makeFilename("chromosome_view.pdf")
   geno.file <- makeFilename("genome_view.pdf")
   cn.file   <- makeFilename("CN_bars.pdf")
   muts.file <- makeFilename("mutations.txt")
   segs.file <- makeFilename("segments.txt")
   robj.extr <- makeFilename("sequenza_extract.RData")
   robj.fit  <- makeFilename("sequenza_cp_table.RData")  
   avg.depth.ratio <- mean(sequenza.extract$gc$adj[, 2])
   assign(x = paste0(sample.id,"_sequenza_extract"), value = sequenza.extract)
   save(list = paste0(sample.id,"_sequenza_extract"), file = robj.extr) 
   if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))){
      stop("Either the cp.table or both cellularity and ploidy argument are required.")
   }
   if (!is.null(cp.table)){
      assign(x = paste0(sample.id,"_sequenza_cp_table"), value = cp.table)
      save(list = paste0(sample.id,"_sequenza_cp_table"), file = robj.fit)       
      cint <- get.ci(cp.table)
      pdf(cp.file)
         cp.plot(cp.table)
         cp.plot.contours(cp.table, add = TRUE, likThresh = c(0.95), col = "red", pch = 20)
         if (!is.null(cellularity) || !is.null(ploidy)) {
            if (is.null(cellularity)) cellularity <- cint$max.cellularity
            if (is.null(ploidy)) ploidy <- cint$max.ploidy
            points(x = ploidy, y = cellularity, pch = 5)
            text(x = ploidy, y = cellularity, 
                 labels = "User selection",
                 pos = 3, offset = 0.5)
         } else {
            cellularity <- cint$max.cellularity
            ploidy <- cint$max.ploidy
         }
      dev.off()
   }
   seg.tab     <- na.exclude(do.call(rbind, sequenza.extract$segments[chromosome.list]))
   mut.tab     <- na.exclude(do.call(rbind, sequenza.extract$mutations[chromosome.list]))
   if (female){
      segs.is.xy <- seg.tab$chromosome == XY["Y"]
      mut.is.xy  <- mut.tab$chromosome == XY["Y"]
   } else{
      segs.is.xy <- seg.tab$chromosome %in% XY
      mut.is.xy  <- mut.tab$chromosome %in% XY
   }

   cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
                            depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio,
                            ratio.priority = ratio.priority, CNn = 2)
   seg.res     <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
   if (!female){
      if (sum(segs.is.xy) >= 1) {
         cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[segs.is.xy], CNt.max = CNt.max,
                               depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                               cellularity = cellularity, ploidy = ploidy,
                               avg.depth.ratio = avg.depth.ratio,
                               ratio.priority = TRUE, CNn = 1)
         seg.xy     <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
         seg.res    <- rbind(seg.res, seg.xy)
      }
   }
   write.table(seg.res, file = segs.file,
               col.names = TRUE, row.names = FALSE, sep = "\t")   

   mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy], CNt.max = CNt.max,
                            depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
                            cellularity = cellularity, ploidy = ploidy,
                            avg.depth.ratio = avg.depth.ratio, CNn = 2)
   mut.res     <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
   if (!female){
      if (sum(mut.is.xy) >= 1) {
         mut.alleles  <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy], CNt.max = CNt.max,
                                  depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                                  cellularity = cellularity, ploidy = ploidy,
                                  avg.depth.ratio = avg.depth.ratio, CNn = 1)
         mut.xy     <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
         mut.res    <- rbind(mut.res, mut.xy)
      }
   }
   write.table(mut.res, file = muts.file,
               col.names = TRUE, row.names = FALSE, sep = "\t")   
   
   pdf(chrw.file)
   for (i in unique(seg.res$chromosome)) {
      
      if (!female && i %in% XY){
         CNn <- 1
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
   pdf(geno.file, height = 5, width = 15)
      genome.view(seg.res)
      genome.view(seg.res, "CN")
   dev.off()
   barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos,
                     CNt = seg.res$CNt)
   cn.sizes <- split(barscn$size, barscn$CNt)
   cn.sizes <- sapply(cn.sizes, 'sum')
   pdf(cn.file)
   barplot(round(cn.sizes / sum(cn.sizes) * 100), names = names(cn.sizes), las = 1,
                      ylab = "Percentage (%)", xlab = "Copy number")
   dev.off()
   
   ## Write down the results.... ploidy etc...
   if (!is.null(cp.table)){
      res.tab <- data.frame(cellularity     = c(cint$confint.cellularity[1], cint$max.cellularity[1], cint$confint.cellularity[2]),
                           ploidy.estimate = c(cint$confint.ploidy[1], cint$max.ploidy[1], cint$confint.ploidy[2]),
                           ploidy.mean.cn  = weighted.mean(x = as.integer(names(cn.sizes)), w = cn.sizes))
      write.table(res.tab, cint.file, col.names = TRUE,
                  row.names = FALSE, sep = "\t")
   }
}
