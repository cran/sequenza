subclonal.matrix <- function(mut.tab, cellularity = seq(0.1, 1, 0.05), ploidy, avg.depth.ratio, mc.cores = 2){
  
  mut.types.list <- lapply(X = 1:nrow(mut.tab),
                           FUN = function(x) { 
                             types.matrix(CNn = mut.tab[x, 'CNn'],
                                          CNt.min = mut.tab[x, 'CNt'],
                                          CNt.max = mut.tab[x, 'CNt'])})
  mut.cloanlity <- function(F, depth.t, types, cellularity, ploidy, avg.depth.ratio) {
    theorethic <- model.points(cellularity = cellularity,
                               ploidy = ploidy,
                               types = types,
                               avg.depth.ratio = avg.depth.ratio)
    max(mufreq.dpois(mufreq = F, mufreq.model = theorethic[, 1], depth.t = depth.t),na.rm = TRUE)
  }
  res <- mclapplyPb (mc.cores = mc.cores, X = 1:nrow(mut.tab),
                     FUN = function (i) {
                             sapply(X = cellularity, FUN = function(x) {
                                mut.cloanlity(F = mut.tab$F[i],
                                           depth.t = mut.tab$good.s.reads[i],
                                           types = mut.types.list[[i]],
                                           cellularity = x, ploidy = ploidy,
                                           avg.depth.ratio = avg.depth.ratio)
                           })
                         })
  res <- do.call(rbind, res)
  res <- res/rowSums(res)
  colnames(res) <- cellularity
  res
}

sequenza2PyClone <- function(mut.tab, seg.cn, sample.id, norm.cn = 2) { 
  mut.tab <- cbind(mut.tab[,c('chromosome', 'n.base', 'good.s.reads','F', 'mutation')], CNt = NA, A = NA, B = NA)
  for (i in 1:nrow(seg.cn)) {
     pos.filt <- mut.tab$chromosome == seg.cn$chromosome[i] & mut.tab$n.base >= seg.cn$start.pos[i] & mut.tab$n.base <= seg.cn$end.pos[i]
     mut.tab[pos.filt, c("CNt", "A", "B")] <- seg.cn[i, c("CNt", "A", "B")]
  }
  id          <- paste(sample.id, mut.tab$chromosome, mut.tab$n.base, sep = ':')
  var.counts  <- round(mut.tab$good.s.reads * mut.tab$F, 0)
  nor.counts  <- mut.tab$good.s.reads - var.counts
  pyclone.tsv <- data.frame(mutation_id = id, ref_counts = nor.counts, var_counts = var.counts,
                            normal_cn = norm.cn,	minor_cn = mut.tab$B, major_cn = mut.tab$A,
                            variant_case = sample.id,	variant_freq = mut.tab$F,	genotype = mut.tab$mutation)
  na.exclude(pyclone.tsv)
}