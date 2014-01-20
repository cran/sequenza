# subclonal.matrix <- function(mut.tab, cellularity = seq(0.1, 1, 0.05), ploidy, avg.depth.ratio, mc.cores = 2){
# 
#   mut.types.list <- lapply(X = 1:nrow(mut.tab),
#                            FUN = function(x) {
#                              types.matrix(CNn = mut.tab[x, 'CNn'],
#                                           CNt.min = mut.tab[x, 'CNt'],
#                                           CNt.max = mut.tab[x, 'CNt'])})
#   mut.cloanlity <- function(F, depth.t, types, cellularity, ploidy, avg.depth.ratio) {
#     theorethic <- model.points(cellularity = cellularity,
#                                ploidy = ploidy,
#                                types = types,
#                                avg.depth.ratio = avg.depth.ratio)
#     max(mufreq.dpois(mufreq = F, mufreq.model = theorethic[, 1], depth.t = depth.t),na.rm = TRUE)
#   }
#   res <- mclapplyPb (mc.cores = mc.cores, X = 1:nrow(mut.tab),
#                      FUN = function (i) {
#                              sapply(X = cellularity, FUN = function(x) {
#                                 mut.cloanlity(F = mut.tab$F[i],
#                                            depth.t = mut.tab$good.s.reads[i],
#                                            types = mut.types.list[[i]],
#                                            cellularity = x, ploidy = ploidy,
#                                            avg.depth.ratio = avg.depth.ratio)
#                            })
#                          })
#   res <- do.call(rbind, res)
#   res <- res/rowSums(res)
#   colnames(res) <- cellularity
#   res
# }

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

VarScan2abfreq <- function(varscan.snp, varscan.copynumber = NULL) {
   
   iupac.nucs     <- setNames(c('A', 'C', 'G', 'GT', 'AC', 'AG', 'CG', 'T', 'AT', 'CT'),
                              c('A', 'C', 'G', 'K', 'M', 'R', 'S', 'T', 'W', 'Y'))
   zygosity.vect  <- setNames(c('hom', 'hom', 'hom', 'het', 'het', 'het', 'het', 'hom', 'het', 'het'),
                              c('A', 'C', 'G', 'K', 'M', 'R', 'S', 'T', 'W', 'Y'))
   varscan.snp  <- varscan.snp[varscan.snp$somatic_status != 'Unknown', ]
   varscan.snp$normal_var_freq <- as.numeric(sub('%', '', varscan.snp$normal_var_freq))/100
   varscan.snp$tumor_var_freq  <- as.numeric(sub('%', '', varscan.snp$tumor_var_freq))/100
   ref.zygosity <- zygosity.vect[varscan.snp$normal_gt]
   AB.germline  <- iupac.nucs[varscan.snp$normal_gt]
   AB.sample    <- rep('.', length(AB.germline))
   
   depth.normal <- varscan.snp$normal_reads1 + varscan.snp$normal_reads2
   depth.sample <- varscan.snp$tumor_reads1 + varscan.snp$tumor_reads1
   depth.ratio  <- depth.sample/depth.normal
   Af <- 1 - varscan.snp$tumor_var_freq
   Bf <- rep(0, length(Af))
   idx <- ref.zygosity == 'het' & Af  < 0.5
   Af[idx] <- 1 - Af[idx]
   idx <- ref.zygosity == 'het'
   Bf[idx] <- 1 - Af[idx]
   idx <- ref.zygosity == 'hom' & varscan.snp$somatic_status == 'Somatic'
   mut <- cbind(as.character(iupac.nucs[varscan.snp$tumor_gt[idx]]),
                as.character(varscan.snp$normal_gt[idx]))
   mut <- sapply(X = 1:sum(idx),
                 FUN = function(x) gsub(x = mut[x, 1],
                                        pattern = mut[x, 2],
                                        replacement = ''))
   mut <- paste0(mut, varscan.snp$tumor_var_freq[idx])
   AB.sample[idx] <- mut
   res <- data.frame(chromosome = as.character(varscan.snp$chrom), n.base = varscan.snp$position,
                     base.ref = as.character(varscan.snp$ref), depth.normal = depth.normal,
                     depth.sample = depth.sample, depth.ratio = depth.ratio,
                     Af = round(Af, 3), Bf = round(Bf, 3), ref.zygosity = ref.zygosity, GC.percent = 50,
                     good.s.reads = round(depth.sample, 2), AB.germline = AB.germline,
                     AB.sample = AB.sample, stringsAsFactors = FALSE)
   normal.pos <- res$ref.zygosity == 'hom' & res$AB.sample == '.'
   res <- res[res$depth.ratio > 0 & !is.infinite(res$depth.ratio) & !normal.pos, ]
   if (!is.null(varscan.copynumber)){
      smart.id <- order(c(1:nrow(varscan.copynumber),1:nrow(varscan.copynumber)+0.5))
      varscan.copynumber$log2_ratio   <- 2^(varscan.copynumber$log2_ratio)
      varscan.copynumber$normal_depth <- round(varscan.copynumber$normal_depth, 0)
      varscan.copynumber$tumor_depth <- round(varscan.copynumber$tumor_depth, 0)
      varscan.copynumber$chrom <- as.character(varscan.copynumber$chrom)
      mat.t <- data.frame(chromosome = c(varscan.copynumber$chrom, varscan.copynumber$chrom)[smart.id],
                          n.base = c(varscan.copynumber$chr_start, varscan.copynumber$chr_stop)[smart.id], base.ref = 'N',
                          depth.normal = c(varscan.copynumber$normal_depth, varscan.copynumber$normal_depth)[smart.id],
                          depth.sample = c(varscan.copynumber$tumor_depth, varscan.copynumber$tumor_depth)[smart.id],
                          depth.ratio  = c(varscan.copynumber$log2_ratio, varscan.copynumber$log2_ratio)[smart.id],
                          Af = 1, Bf = 0, ref.zygosity = 'hom',
                          GC.percent = c(varscan.copynumber$gc_content, varscan.copynumber$gc_content)[smart.id],
                          stringsAsFactors = FALSE)
      mat.t <- cbind(mat.t, good.s.reads = mat.t$depth.sample, AB.germline = 'N', AB.sample = '.')
      chrom.order <- unique(mat.t$chromosome)
      l.cnv <- split(mat.t, mat.t$chromosome)
      l.snp <- split(res, res$chromosome)     
      for (i in names(l.snp)){
         tab.i <- rbind(l.cnv[[i]], l.snp[[i]])
         tab.i <- tab.i[order(tab.i$n.base), ]
         #idx.i <- diff(tab.i$n.base)
         #dups.i  <- which(idx.i == 0)
         #snp.i <- which(tab.i$GC.percent == -100)
         l.cnv[[i]] <- tab.i
      }
      res <- do.call(rbind, l.cnv[chrom.order])
   }
   res
}