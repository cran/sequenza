read.seqz <- function(file, n_lines = NULL, col_types = "ciciidddcddccc",
    chr_name = NULL, buffer = 33554432, parallel = 1,
    col_names = c("chromosome", "position", "base.ref", "depth.normal",
                  "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
                  "GC.percent", "good.reads", "AB.normal", "AB.tumor",
                  "tumor.strand"), ...) {

    if (is.null(n_lines)) {
        skip <- 1
        n_max <- Inf
    } else {
        n_lines <- round(sort(n_lines), 0)
        skip <- n_lines[1]

        n_max <- n_lines[2] - skip + 1
    }
    if (!is.null(chr_name)) {
        chr_name <- as.character(chr_name)
        tbi <- file.exists(paste(file, "tbi", sep = "."))
        if (tbi) {
            read.seqz.tbi(file, split_chr_coord(chr_name),
                col_names, col_types)
        } else {
            read.seqz.chr(file, chr_name = chr_name, col_types = col_types,
                col_names = col_names, buffer = buffer, parallel = parallel)
        }
    } else {
        read_tsv(file = file, col_types = col_types, skip = skip,
            n_max = n_max, col_names = col_names, progress = FALSE, ...)
    }
}

read.seqz.chr <- function(file, chr_name, col_names,
    col_types, buffer, parallel) {
    con <- gzfile(file, "rb")
    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, chr_name, col_names, col_types) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, skip = 0, n_max = Inf,
            col_names = col_names, progress = FALSE)
        x[x$chromosome == chr_name, ]
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, chr_name = chr_name,
        col_names = col_names, col_types = col_types, CH.MAX.SIZE = buffer,
        parallel = parallel)
    close(con)
    res
}

read.seqz.tbi <- function(file, chr_name, col_names, col_types) {
    res <- tabix.read(file, chr_name)
    res <- read_tsv(file = paste(mstrsplit(res), collapse = "\n"),
        col_types = col_types, skip = 0, n_max = Inf,
        col_names = col_names, progress = FALSE)
}
