mut.fractions <- function(AB.tumor, Af, tumor.strand) {
    F = 1 - Af
    base.mut <- lapply(X = AB.tumor, FUN = function(x) {
        unlist(strsplit(as.character(x), split = "[:]"))
    })
    base.fw  <- lapply(X = tumor.strand, FUN = function(x) {
        unlist(strsplit(as.character(x), split = "[:]"))
    })
    frequencify <- function (x) {
        base.name <- substr(unlist(x), 1, 1)
        base.val  <- as.numeric(substr(unlist(x), 2, nchar(x)))
        setNames(base.val, base.name)
    }
    base.freqs <- lapply(X = base.mut, FUN = frequencify)
    fw.freqs   <- lapply(X = base.fw,  FUN = frequencify)
    n.base.mut <- do.call(c, lapply(X = base.mut, FUN = length))
    max.fq <- function (x) {
        freq.rel <- base.freqs[[x]] / F[x]
        f.max    <- which.max(freq.rel)
        c(freq.rel[f.max], names(base.freqs[[x]])[f.max],
            base.freqs[[x]][f.max], fw.freqs[[x]][f.max])
    }
    max.freqs  <- do.call(rbind, lapply(1:length(F), max.fq))
    data.frame(base.count = as.integer(n.base.mut),
        maj.base.freq = as.numeric(max.freqs[, 1]),
        base = as.character(max.freqs[, 2]),
        freq = as.numeric(max.freqs[, 3]),
        fw.freq = as.numeric(max.freqs[, 4]))
}

mutation.table <- function(seqz.tab, mufreq.treshold = 0.15,
    min.reads = 40, min.reads.normal = 10, max.mut.types = 3,
    min.type.freq = 0.9, min.fw.freq = 0, segments = NULL) {
    chroms <- unique(seqz.tab$chromosome)
    hom.filt <- seqz.tab$zygosity.normal == "hom" &
       seqz.tab$AB.tumor != "."
    seqz.tab <- seqz.tab[hom.filt, ]
    reads.filt <- seqz.tab$good.reads >= min.reads &
        seqz.tab$depth.normal >= min.reads.normal
    seqz.tab <- seqz.tab[reads.filt, ]
    mufreq.filt <- seqz.tab$Af <= (1 - mufreq.treshold)
    seqz.tab <- seqz.tab[mufreq.filt, ]
    if (!is.null(segments)) {
        for (i in 1:nrow(segments)) {
            pos.filt <- seqz.tab$chromosome == segments$chromosome[i] &
                seqz.tab$position >= segments$start.pos[i] &
                seqz.tab$position <= segments$end.pos[i]
            seqz.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
        }
    }
    seqz.dummy   <- data.frame(chromosome = chroms, position = 1,
        GC.percent = NA, good.reads = NA, adjusted.ratio = NA,
        F = 0, mutation = "NA", stringsAsFactors = FALSE)
    if (nrow(seqz.tab) >= 1) {
        mu.fracts   <- mut.fractions(AB.tumor = seqz.tab$AB.tumor,
            Af = seqz.tab$Af, tumor.strand = seqz.tab$tumor.strand)
        mufreq.filt <- mu.fracts$freq >= mufreq.treshold
        type.filt   <- mu.fracts$base.count <= max.mut.types
        prop.filt   <- mu.fracts$maj.base.freq >= min.type.freq
        if (!is.na(min.fw.freq)) {
            fw.2 <- 1 - min.fw.freq
            fw.2 <- sort(c(fw.2, min.fw.freq))
            fw.filt     <- mu.fracts$fw.freq > fw.2[1] &
                mu.fracts$fw.freq < fw.2[2]
            mufreq.filt <- mufreq.filt & type.filt  & prop.filt & fw.filt
        } else {
            mufreq.filt <- mufreq.filt & type.filt  & prop.filt
        }
        mut.type    <- paste(seqz.tab$AB.normal, mu.fracts$base, sep = ">")
        seqz.tab     <- seqz.tab[, c("chromosome", "position", "GC.percent",
            "good.reads", "adjusted.ratio")]
        seqz.tab     <- cbind(seqz.tab, F = mu.fracts$freq,
            mutation = mut.type)
        rbind(seqz.tab[mufreq.filt, ], seqz.dummy)
    } else {
        seqz.dummy
   }
}
