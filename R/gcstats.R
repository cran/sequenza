gc.sample.stats <- function(file, col_types = "c--dd----d----",
    buffer = 33554432, parallel = 2L, verbose = TRUE) {

    con <- gzfile(file, "rb")

    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, col_types) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, col_names = FALSE,
            skip = 0, n_max = Inf, progress = FALSE)
        u_chr <- unique(x[, 1])
        n_chr <- table(x[, 1])
        gc1 <- lapply(split(x[, 2], x[, 4]), table)
        gc2 <- lapply(split(x[, 3], x[, 4]), table)
        if (verbose){
            message(".", appendLF = FALSE)
        }
        list(unique = u_chr, lines = n_chr, gc_nor = gc1, gc_tum = gc2)
    }
    if (verbose){
        message("Collecting GC information ", appendLF = FALSE)
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, col_types = col_types,
        CH.MAX.SIZE = buffer, parallel = parallel)
    close(con)
    if (verbose){
        message(" done\n")
    }
    unfold_gc(res, stats = TRUE)
}

unfold_gc <- function(x, stats = TRUE) {
    gc_norm <- get_gc(x[, "gc_nor"])
    gc_tum <- get_gc(x[, "gc_tum"])
    if (stats) {
        ord_chrom <- unique(Reduce("c", Reduce("c", x[, "unique"])))
        stats_chrom <- Reduce("c", x[, "lines"])
        stats_chrom <- sapply(splash_table(x[, "lines"]), sum)
        stats_chrom <- stats_chrom[ord_chrom]
        stats_start <- cumsum(c(1, stats_chrom[-length(stats_chrom)]))
        stats_end   <- stats_start + stats_chrom - 1
        stats_chrom <- data.frame(chr = ord_chrom, n_lines = stats_chrom,
            start = stats_start, end = stats_end)

        list(file.metrics = stats_chrom, normal = gc_norm, tumor = gc_tum)
    } else {
        list(normal = gc_norm, tumor = gc_tum)
    }
}

splash_table <- function(lis_obj){
    lis_obj <- Reduce("c", lis_obj)
    split(lis_obj, names(lis_obj))
}

get_gc <- function(gc_col) {
    sort_char <- function(x) {
        as.character(sort(as.numeric(x)))
    }
    all_depths <- splash_table(gc_col)
    all_depths <- lapply(all_depths, FUN = function(x) {
        sapply(splash_table(x), sum)
    })
    names_gc <- sort_char(names(all_depths))
    all_depths <- all_depths[names_gc]
    names_depths <- sort_char(unique(Reduce("c", lapply(all_depths, names))))
    n <- do.call(rbind, lapply(all_depths, FUN = function(x, names_depths) {
            res <- x[names_depths]
            names(res) <- names_depths
            res
        },
        names_depths = names_depths))
    n[is.na(n)] <- 0
    list(gc = as.numeric(names_gc), depth = as.numeric(names_depths), n = n)
}

median_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.median(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

mean_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.mean(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

depths_gc <- function(depth_n, depth_t, gc) {
    gc_nor <- lapply(split(depth_n, gc), table)
    gc_tum <- lapply(split(depth_t, gc), table)
    list(gc_nor = gc_nor, gc_tum = gc_tum)
}
