read.abfreq <- function (file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
    colClasses = c('character', 'integer', 'character', 'integer', 
      'integer', 'numeric', 'numeric', 'numeric', 'character', 
      'numeric', 'numeric', "character", "character"), chr.name = NULL, n.lines = NULL, ...) {
   if (!is.null(n.lines) & is.null(chr.name)) fast <-  FALSE
   if(fast && nrows == -1) {
    if(gz) {
       if (!is.null(chr.name)) {
          wc <- system(paste(paste('gzip -d -c | grep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('gzip -d -c', file, '| wc'), intern = TRUE)
       }
    } else {
       if (!is.null(chr.name)) {
          wc <- system(paste(paste('grep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('wc', file), intern = TRUE)
       }
    }
    if (is.null(chr.name)) {
       wc <- sub("^ +", "", wc)
       wc <- strsplit(wc, ' ')[[1]][1]
    }
    nrows <- max(as.integer(wc), 1)
    message('Reading ', nrows, ' lines...')
  }
   if (!is.null(chr.name)) {
      if (gz) {
         grep.part <- paste("gzip -d -c | grep '^", chr.name, "\t'", sep = "")
      } else {
         grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
      }
      abf.data   <- read.delim(pipe(paste(grep.part, file, sep = " ")), nrows = nrows, colClasses = colClasses, header = FALSE, ...)
      if (header == TRUE) {
         head       <- colnames(read.table(file, header = TRUE, nrows = 1 ))
         colnames(abf.data) <- head
      }
      abf.data
   } else {
      if (!is.null(n.lines)){
         if (!is.numeric(n.lines) | length(n.lines) != 2) stop("n.lines must be a vector of 2 integers")
         n.lines <- round(sort(n.lines), 0)
         if (header == TRUE) {
            n.lines <- n.lines + 1
         }
         if(gz) {
            abf.data <- read.delim(pipe(paste("gzip -d -c", file,"| sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'")),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE,...)
         }  else{
            abf.data <- read.delim(pipe(paste("sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'", file)),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE, ...)
         }
         if (header == TRUE) {
            head  <- colnames(read.table(file, header = TRUE, nrows = 1 ))
            colnames(abf.data) <- head
         }
         abf.data
      } else {
         read.delim(file, nrows = nrows, colClasses = colClasses, header = header, ...)
      }
   }
}

read.acgt <- function (file, colClasses = c('character', 'integer', 'character', 'integer', 
                                            'integer', 'integer', 'integer', 'integer'), ...) {
   read.abfreq(file = file , colClasses = colClasses, ...)
}

gc.norm <- function (x, gc) {
   dr.by.gc <- split(x, gc)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.median <- sapply(dr.by.gc, median, na.rm = TRUE)
   dr.by.gc.mean <- sapply(dr.by.gc, mean, na.rm = TRUE)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)),
        raw.mean = dr.by.gc.mean, raw.median = dr.by.gc.median)
}

gc.sample.stats <- function (file, gz = TRUE) {
   colClasses = c('character', 'numeric', 'numeric')
   if (gz) {
      abf.data <- read.delim(pipe(paste('gzip -d -c', file, '| cut -f 1,6,10')), colClasses = colClasses)
   } else {
      abf.data <- read.delim(pipe(paste('cut -f 1,6,10', file)), colClasses = colClasses)
   }
   gc.stats <- gc.norm(x = abf.data$depth.ratio,
                       gc = abf.data$GC.percent)
   chr.ord  <- unique(abf.data$chromosome)
   chr.dim  <- lapply(X = split(abf.data$chromosome, abf.data$chromosome), FUN = length)
   chr.dim  <- data.frame(chr = chr.ord, n.lines = do.call(rbind,chr.dim[chr.ord]))
   chr.dim$start <- cumsum(c(1, chr.dim$n.lines[-length(chr.dim$n.lines)]))
   chr.dim$end   <- chr.dim$start + chr.dim$n.lines - 1
   gc.stats$file.metrics <- chr.dim
   gc.stats
}

windowValues <- function(x, positions, chromosomes, window = 1e6, overlap = 0, verbose = TRUE,
                          weight = rep.int( x = 1, times = length(x)), start.coord = 1) {
   results  <- list()
   weight   <- sqrt(weight)
   xw       <- x * weight
   overlap  <- as.integer(overlap)
   window.o <- window - round(window * (overlap / (overlap + 1)), 0)
   chr.order <- unique(chromosomes)
   mat.w    <- data.frame(chr = chromosomes, pos = positions, x = x, xw = x * weight,
                          w = weight, stringsAsFactors = FALSE)
   mat.w    <- split(mat.w, mat.w$chr)
   do.windows <- function(x.i, w.i, xw.i, breaks, overlap){
      coords    <- data.frame(start = breaks[-length(breaks)], end = breaks[-1])
      if (overlap > 0) {
      start.w   <- 1:nrow(coords)
      end.w     <- start.w + overlap
      end.w[end.w > max(start.w)] <- max(start.w)
      coords    <- data.frame(start = coords[start.w, 1], end = coords[end.w, 2])
      idx.merge <- apply(cbind(start.w, end.w), 1, unique)
      x.i       <- lapply(1:nrow(coords), function(a) do.call(c, x.i[idx.merge[[a]]]))
      w.i       <- lapply(1:nrow(coords), function(a) do.call(c, w.i[idx.merge[[a]]]))
      xw.i      <- lapply(1:nrow(coords), function(a) do.call(c, xw.i[idx.merge[[a]]]))
      }
      quartiles <- do.call(rbind, 
                           lapply(X = x.i, FUN = function(a) quantile(a, probs = c(0.25, 0.75),
                                                                      na.rm = TRUE)))
      sum.w     <- sapply(X = w.i, FUN = function(a) sum(a, na.rm = TRUE))
      sum.xw    <- sapply(X = xw.i, FUN = function(a) sum(a, na.rm = TRUE))
      size      <- sapply(X = x.i, FUN = length)
      data.frame(coords, mean = sum.xw/sum.w, q0 = quartiles[,1],
                 q1 = quartiles[,2], N = size, row.names = 1:length(size))            
   }
   for (i in 1:length(mat.w)) {
      range.pos            <- range(mat.w[[i]]$pos, na.rm = TRUE)
      if (!is.null(start.coord)) {
         range.pos[1] <- as.integer(start.coord)
      }
      beam.coords          <- seq(range.pos[1], range.pos[2], by = window.o)
      if (max(beam.coords) != range.pos[2] ) {
         beam.coords <- c(beam.coords, range.pos[2])
      }
      f.windows <- cut(x = mat.w[[i]]$pos, breaks = beam.coords)
      xw      <- split(x = mat.w[[i]]$xw, f = f.windows)
      w       <- split(x = mat.w[[i]]$w, f = f.windows)
      x       <- split(x = mat.w[[i]]$x, f = f.windows)
      if (overlap > 0 ) {
         if (verbose) {
            cat(paste("chromosome:", names(mat.w)[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window, "overlapping windows:", overlap, "\n",sep=" "))
         }
      } else {
         if (verbose) {
            cat(paste("chromosome:", names(mat.w)[i], "from:", range.pos[1],"to:", range.pos[2], "window:", window,"\n",sep=" "))
         }
      }
      
      results[[i]] <- do.windows(x.i = x, w.i = w, xw.i = xw,
                                 breaks = beam.coords, overlap = overlap)
   }
   results <- results[as.factor(chr.order)]
   names(results) <- names(chr.order)
   results
}

get.ci <- function(cp.table, interval = 0.95) {
   results  <- list()
   maxs.x   <- apply(cp.table$z, 1, FUN = max)
   maxs.y   <- apply(cp.table$z, 2, FUN = max)
   max.xy   <- which(cp.table$z == max(cp.table$z), arr.ind = TRUE)
   val.95.x <- quantile(maxs.x, prob = interval, na.rm = TRUE)
   val.95.y <- quantile(maxs.y, prob = interval, na.rm = TRUE)
   values.x <- cp.table$x[maxs.x >= val.95.x]
   values.y <- cp.table$y[maxs.y >= val.95.y]
   up.x     <- max(values.x)
   low.x    <- min(values.x) 
   max.x    <- cp.table$x[which.max(maxs.x)]
   up.y     <- max(values.y)
   low.y    <- min(values.y) 
   max.y    <- cp.table$y[which.max(maxs.y)]
   results$values.x  <- cbind(x = cp.table$x, y = maxs.x)
   results$confint.x <- c(low.x, up.x)
   results$max.x     <- max.x
   results$values.y  <- cbind(x = maxs.y, y = cp.table$y)
   results$confint.y <- c(low.y, up.y)
   results$max.y     <- max.y
   results
}

# merge.baf.ratio <- function(baf.segments, ratio.segments) {
#    baf.table    <- lapply(1:length(baf.segments), FUN = function(x) cbind(chromosome = names(baf.segments)[x], 
#                                                                         as.data.frame(do.call(rbind, baf.segments[[x]]))))
#    ratio.table  <- lapply(1:length(ratio.segments), FUN = function(x) cbind(chromosome = names(ratio.segments)[x],
#                                                                           as.data.frame(do.call(rbind, ratio.segments[[x]]))))
#    baf.table    <- do.call(rbind, baf.table)
#    ratio.table  <- do.call(rbind, ratio.table)
#    baf.table$chromosome <- as.character(baf.table$chromosome)
#    ratio.table$chromosome <- as.character(ratio.table$chromosome)
#    baf.index    <- sapply(1:nrow(baf.table), FUN = function(x) paste(baf.table[x, 1:2], collapse ="_"))
#    ratio.index  <- sapply(1:nrow(ratio.table), FUN = function(x) paste(ratio.table[x, 1:2], collapse = "_"))
#    baf.table    <- cbind(index = baf.index, baf.table[, c(5,7)])
#    ratio.table  <- cbind(index = ratio.index, ratio.table[, c(1,2,3,5,7)])
#    colnames(ratio.table) <- c("index", "chromosome", "start", "end", "ratio", "N.depth")
#    colnames(baf.table) <- c("index", "Bf", "N.BAF")
#    merged.table <- merge(ratio.table, baf.table, all = TRUE, sort = FALSE)
#    merged.table <- merged.table[, -1]

#    merged.table
# }


mut.fractions <- function(AB.sample, Af) {
  F = 1 - Af
   base.mut <- lapply(X = AB.sample, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   frequencify <- function (x) {
      base.name <- substr(unlist(x), 1, 1)
      base.val  <- as.numeric(substr(unlist(x), 2, nchar(x)))
      setNames(base.val, base.name)
   }
   base.freqs <- lapply(X = base.mut, FUN = frequencify)   
   n.base.mut <- do.call(c, lapply(X = base.mut, FUN = length))
   max.fq <- function (x) {
      freq.rel <- base.freqs[[x]] / F[x]
      f.max    <- which.max(freq.rel)
      c(freq.rel[f.max], names(base.freqs[[x]])[f.max], base.freqs[[x]][f.max])
   }
   max.freqs  <- do.call(rbind, lapply(1:length(F), max.fq))
   data.frame(base.count = as.integer(n.base.mut), maj.base.freq = as.numeric(max.freqs[, 1]),
              base = as.character(max.freqs[,2]), freq = as.numeric(max.freqs[,3]))
}

mutation.table <- function(abf.tab, mufreq.treshold = 0.15, min.reads = 40, max.mut.types = 3,
                           min.type.freq = 0.9, segments = NULL) {
   chroms      <- unique(abf.tab$chromosome)
   hom.filt    <- abf.tab$ref.zygosity == 'hom'
   abf.tab     <- abf.tab[hom.filt, ]
   reads.filt  <- abf.tab$good.s.reads >= min.reads
   abf.tab     <- abf.tab[reads.filt, ]
   mufreq.filt <- abf.tab$Af <= (1 - mufreq.treshold)
   abf.tab     <- abf.tab[mufreq.filt, ]
   if (!is.null(segments)) {
      for (i in 1:nrow(segments)) {
         pos.filt <- abf.tab$chromosome == segments$chrom[i] & abf.tab$n.base >= segments$start.pos[i] & abf.tab$n.base <= segments$end.pos[i]
         abf.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
      }
   }
   abf.dummy   <- data.frame(chromosome = chroms, n.base = 1, GC.percent = NA, good.s.reads = NA,
                             adjusted.ratio = NA, F = 0, mutation = 'NA', stringsAsFactors= FALSE)   
   if (nrow(abf.tab) >= 1) {
      mu.fracts   <- mut.fractions(AB.sample = abf.tab$AB.sample, Af = abf.tab$Af)
      mufreq.filt <- mu.fracts$freq >= mufreq.treshold
      type.filt   <- mu.fracts$base.count <= max.mut.types
      prop.filt   <- mu.fracts$maj.base.freq <= min.type.freq
      mut.type    <- paste(abf.tab$AB.germline, mu.fracts$base, sep = '>')   
      abf.tab     <- abf.tab[,c('chromosome', 'n.base', 'GC.percent', 'good.s.reads', 'adjusted.ratio')]
      abf.tab     <- cbind(abf.tab, F = mu.fracts$freq, mutation = mut.type)
      rbind(abf.tab[mufreq.filt, ], abf.dummy)
   } else {
      abf.dummy
   }
}

find.breaks <- function(abf.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), verbose = FALSE, ...) {
   chromosome <- gsub(x = abf.baf$chromosome, pattern = "chr", replacement = "")
   logR = data.frame(chrom = chromosome, 
                     pos = abf.baf$n.base,
                     s1 = log2(abf.baf$adjusted.ratio))
   BAF = data.frame(chrom = chromosome, 
                    pos = abf.baf$n.base,
                    s1 = abf.baf$Bf)
   logR.wins <- copynumber::winsorize(logR, verbose = verbose)
   allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF, baf.thres = baf.thres,
                       verbose = verbose, gamma = gamma, kmin = kmin, ...)
    if (length(grep("chr", abf.baf$chromosome)) > 0) { 
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
    allele.seg[allele.seg$end.pos - allele.seg$start.pos != 0,
               c("chrom", "start.pos", "end.pos")]
}

segment.breaks <- function(abf.tab, breaks) {
   w.r     <- sqrt(abf.tab$depth.sample)
   rw      <- abf.tab$adjusted.ratio * w.r
   w.b     <- sqrt(abf.tab$good.s.reads)
   bw      <- abf.tab$Bf * w.b
   abf.tab <- cbind(abf.tab[, c("chromosome", "n.base", "ref.zygosity")],
                    rw = rw, w.r = w.r, bw = bw, w.b = w.b)
   chr.order <- unique(abf.tab$chromosome)
   abf.tab <- split(abf.tab, f = abf.tab$chromosome)
   segments <- list()
   for (i in 1:length(abf.tab)) {
      abf.b.i     <- abf.tab[[i]][abf.tab[[i]]$ref.zygosity == 'het', ]
      breaks.i    <- breaks[breaks$chrom == names(abf.tab)[i], ]
      nb          <- nrow(breaks.i)
      breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,c("start.pos", "end.pos")], f = 1:nb))
      fact.r.i    <- cut(abf.tab[[i]]$n.base, breaks.vect)
      fact.b.i    <- cut(abf.b.i$n.base, breaks.vect)
      seg.i.s.r   <- sapply(X = split(abf.tab[[i]]$w.r, f = fact.r.i), FUN = length)
      seg.i.s.b   <- sapply(X = split(abf.b.i$w.b, f = fact.b.i), FUN = length)      
      seg.i.rw    <- sapply(X = split(abf.tab[[i]]$rw, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
      seg.i.w.r   <- sapply(X = split(abf.tab[[i]]$w.r, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
      seg.i.bw    <- sapply(X = split(abf.b.i$bw, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
      seg.i.w.b   <- sapply(X = split(abf.b.i$w.b, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
      segments.i <- data.frame(chromosome  = names(abf.tab)[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                               end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.bw/seg.i.w.b, N.BAF = seg.i.s.b,
                               depth.ratio = seg.i.rw/seg.i.w.r, N.ratio = seg.i.s.r, stringsAsFactors = FALSE)
      segments[[i]] <- segments.i[seq(from = 1, to = nrow(segments.i), by = 2),]
   }
   segments <- do.call(rbind, segments[as.factor(chr.order)])
   row.names(segments) <- 1:nrow(segments)
   segments
}
