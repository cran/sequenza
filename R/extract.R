sequenza.extract <- function(file, window = 1e6, overlap = 1,
    gamma = 80, kmin = 10, gamma.pcf = 140, kmin.pcf = 40,
    mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
    breaks = NULL, breaks.method = "het", assembly = "hg19",
    weighted.mean = TRUE, normalization.method = "mean",
    ignore.normal = FALSE, parallel = 1, gc.stats = NULL,
    segments.samples = FALSE){

    if (is.null(gc.stats)) {
        gc.stats <- gc.sample.stats(file, verbose = verbose,
            parallel = parallel)
    }
    if (normalization.method == "mean") {
        gc.normal.vect <- mean_gc(gc.stats$normal)
        gc.tumor.vect  <- mean_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.mean(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.mean(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    } else {
        gc.normal.vect <- median_gc(gc.stats$normal)
        gc.tumor.vect  <- median_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.median(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.median(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    }
    windows.baf   <- list()
    windows.ratio <- list()
    windows.raw_ratio <- list()
    windows.normal <- list()
    windows.tumor <- list()
    windows.n_normal <- list()
    windows.n_tumor <- list()
    mutation.list <- list()
    segments.list <- list()
    segments_samples.list <- list()
    norm.gc.list <- list()
    if (is.null(dim(breaks))) {
        breaks <- NULL
    }
    chr.vect <- as.character(gc.stats$file.metrics$chr)
    if (is.null(chromosome.list)) {
        chromosome.list <- chr.vect
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }
    for (chr in chromosome.list) {
        if (verbose) {
            message("Processing ", chr, ":", appendLF = TRUE)
        }
        tbi <- file.exists(paste0(file, ".tbi"))
        if (tbi) {
            seqz.data   <- read.seqz(file, chr_name = chr)
        } else {
            file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
            seqz.data   <- read.seqz(file, n_lines = c(file.lines$start,
                file.lines$end))
        }

        norm_tumor_depth <- seqz.data$depth.tumor /
            gc.tumor.vect[as.character(seqz.data$GC.percent)]
        norm_normal_depth <- seqz.data$depth.normal /
            gc.normal.vect[as.character(seqz.data$GC.percent)]
        norm.gc.stats <- depths_gc(
            depth_n = round(norm_normal_depth * avg_nor_depth, 0),
            depth_t = round(norm_tumor_depth * avg_tum_depth, 0),
            gc = seqz.data$GC.percent)
        if (ignore.normal) {
            seqz.data$adjusted.ratio <- round(norm_tumor_depth, 3)
        } else {
            seqz.data$adjusted.ratio <- round(
                norm_tumor_depth / norm_normal_depth, 3)
        }
        if (segments.samples == TRUE) {
            breaks_normal_chr <- breaks_full(
                data = data.frame(chromosome = seqz.data$chromosome,
                                  position = seqz.data$position,
                                  adjusted.ratio = norm_normal_depth,
                                  singsAsFactors = FALSE),
                gamma = gamma.pcf, kmin = kmin.pcf, assembly = assembly,
                breaks.het = NULL)
            breaks_tumor_chr <- breaks_full(
               data = data.frame(chromosome = seqz.data$chromosome,
                                 position = seqz.data$position,
                                 adjusted.ratio = norm_tumor_depth,
                                 singsAsFactors = FALSE),
               gamma = gamma.pcf, kmin = kmin.pcf, assembly = assembly,
               breaks.het = NULL)
        } else {
           breaks_normal_chr <- NULL
           breaks_tumor_chr <- NULL
        }
        seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap,
            weight = seqz.data$depth.normal)
        seqz.n.win <- windowValues(x = seqz.data$depth.normal / avg_nor_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.t.win <- windowValues(x = seqz.data$depth.tumor / avg_tum_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.r_r.win <- windowValues(x = seqz.data$depth.ratio / (
                avg_tum_depth / avg_nor_depth),
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap,
            weight = seqz.data$depth.normal)
        seqz.n_n.win <- windowValues(x = norm_normal_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.n_t.win <- windowValues(x = norm_tumor_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)


        seqz.hom <- seqz.data$zygosity.normal == "hom"
        seqz.het <- seqz.data[!seqz.hom, ]
        het.filt <- seqz.het$good.reads >= min.reads.baf
        seqz.het <- seqz.het[het.filt, ]
        het_ok <- nrow(seqz.het) > 0
        if (is.null(breaks)) {
            breaks_chr <- NULL
        } else {
            breaks_chr <- breaks[breaks$chrom == chr, ]
        }
        if (het_ok) {
            seqz.b.win <- windowBf(Af = seqz.het$Af, Bf = seqz.het$Bf,
                good.reads = seqz.het$good.reads,
                chromosomes = seqz.het$chromosome,
                positions = seqz.het$position, conf = 0.95,
                window = window, overlap = overlap)
        } else {
            seqz.b.win <- list()
            seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position,
                na.rm = TRUE), end = max(seqz.data$position, na.rm = TRUE),
                mean = 0, q0 = 0,  q1 = 0, N = 1)
        }
        if (het_ok) {
            breaks_chr <- extract_breaks(data = seqz.data, data_het = seqz.het,
                ratio = seqz.r.win, baf = seqz.b.win,
                gamma = gamma, kmin = kmin, breaks = breaks_chr,
                gamma.pcf = gamma.pcf, kmin.pcf = kmin.pcf,
                assembly = assembly, chromosome = chr,
                method = breaks.method)
        } else {
            if (breaks.method == "full") {
                breaks_chr <- extract_breaks(data = seqz.data,
                    data_het = seqz.het, ratio = seqz.r.win, baf = seqz.b.win,
                    gamma = gamma, kmin = kmin,
                    gamma.pcf = gamma.pcf, kmin.pcf = kmin.pcf,
                    assembly = assembly, chromosome = chr,
                    method = breaks.method)
            }
        }
        if (class(breaks_chr) == "try-error") {
           breaks_chr <- NULL
        }
        if (is.null(breaks_chr) || nrow(breaks_chr) == 0 ||
            length(breaks_chr) == 0) {
            breaks_chr <- data.frame(chrom = chr,
                start.pos = min(seqz.data$position, na.rm = TRUE),
                end.pos = max(seqz.data$position, na.rm = TRUE))
        }
        seg.s1 <- segment.breaks(seqz.tab = seqz.data, breaks = breaks_chr,
            min.reads.baf = min.reads.baf, weighted.mean = weighted.mean)

        mut.tab   <- mutation.table(seqz.data,
            mufreq.treshold = mufreq.treshold,
            min.reads = min.reads, min.reads.normal = min.reads.normal,
            max.mut.types = max.mut.types, min.type.freq = min.type.freq,
            min.fw.freq = min.fw.freq, segments = seg.s1)

        windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
        windows.raw_ratio[[which(chromosome.list == chr)]] <- seqz.r_r.win[[1]]
        windows.normal[[which(chromosome.list == chr)]] <- seqz.n.win[[1]]
        windows.tumor[[which(chromosome.list == chr)]] <- seqz.t.win[[1]]
        windows.n_normal[[which(chromosome.list == chr)]] <- seqz.n_n.win[[1]]
        windows.n_tumor[[which(chromosome.list == chr)]] <- seqz.n_t.win[[1]]
        windows.baf[[which(chromosome.list == chr)]]   <- seqz.b.win[[1]]
        segments.list[[which(chromosome.list == chr)]] <- seg.s1
        mutation.list[[which(chromosome.list == chr)]] <- mut.tab
        norm.gc.list[[which(chromosome.list == chr)]] <- norm.gc.stats
        segments_samples.list[[which(chromosome.list == chr)]] <- list(
            normal = breaks_normal_chr, tumor = breaks_tumor_chr)

        if (verbose) {
            message('   ', nrow(mut.tab), ' variant calls.', appendLF = TRUE)
            message('   ', nrow(seg.s1), ' copy-number segments.', appendLF = TRUE)
            message('   ', nrow(seqz.het), ' heterozygous positions.', appendLF = TRUE)
            message('   ', sum(seqz.hom), ' homozygous positions.', appendLF = TRUE)
        }
    }
    names(windows.baf)   <- chromosome.list
    names(windows.ratio) <- chromosome.list
    names(windows.raw_ratio) <- chromosome.list
    names(windows.normal) <- chromosome.list
    names(windows.tumor) <- chromosome.list
    names(windows.n_normal) <- chromosome.list
    names(windows.n_tumor) <- chromosome.list
    names(mutation.list) <- chromosome.list
    names(segments.list) <- chromosome.list
    names(segments_samples.list) <- chromosome.list


    gc_norm <- unfold_gc(do.call(rbind, norm.gc.list), stats = FALSE)

    if (normalization.method == "mean") {
        avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth,
            w = colSums(gc_norm$tumor$n))
        avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth,
            w = colSums(gc_norm$normal$n))
    } else {
        avg_tum_ndepth <- weighted.median(x = gc_norm$tumor$depth,
            w = colSums(gc_norm$tumor$n))
        avg_nor_ndepth <- weighted.median(x = gc_norm$normal$depth,
            w = colSums(gc_norm$normal$n))
    }
    if (ignore.normal) {
        avg_depth_ratio <- avg_tum_ndepth / avg_tum_depth
    } else {
        avg_depth_ratio <- (avg_tum_ndepth / avg_tum_depth) /
            (avg_nor_ndepth / avg_nor_depth)
    }

    list(BAF = windows.baf, ratio = windows.ratio,
        raw_ratio = windows.raw_ratio,
        depths =  list(
            raw = list(normal = windows.normal, tumor = windows.tumor),
            norm = list(normal = windows.n_normal, tumor = windows.n_tumor)),
        mutations = mutation.list, segments = segments.list,
        chromosomes = chromosome.list, gc = gc.stats,
        gc_norm = gc_norm, avg.depth.ratio = avg_depth_ratio,
        avg.depth.tumor = avg_tum_depth, avg.depth.normal = avg_nor_depth,
        segments_samples = segments_samples.list)
}
