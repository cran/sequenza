### R code from vignette source 'sequenzaSNP.Rnw'

###################################################
### code chunk number 1: loadingSequenza
###################################################
library(sequenza)


###################################################
### code chunk number 2: loadingFromCopynumber
###################################################
data(BAF)
data(logR)


###################################################
### code chunk number 3: SNPdataPrep
###################################################
sample.i <- data.frame(chromosome = BAF$chrs, position = BAF$pos, 
                      Bf = BAF$S1, adjusted.ratio = logR$S1,
                      depth.tumor = 1, good.reads = 1,
                      zygosity.normal = 'hom', stringsAsFactors = FALSE)


###################################################
### code chunk number 4: logRMeansub
###################################################

sample.i$adjusted.ratio <- 2^(sample.i$adjusted.ratio/0.55)



###################################################
### code chunk number 5: hetFind
###################################################
het.lim <- 0.2
is.het <- sample.i$Bf >= het.lim & sample.i$Bf <= 1 - het.lim
sample.i$zygosity.normal[is.het] <- 'het'
sample.i$Bf[sample.i$Bf >= 0.5] <- 1 - sample.i$Bf[sample.i$Bf >= 0.5]
sample.het.i <- sample.i[is.het, ]


###################################################
### code chunk number 6: logRWin
###################################################
 
snp.r.win <- windowValues(x = sample.i$adjusted.ratio,
                          positions = sample.i$position,
                          chromosomes = sample.i$chromosome,
                          window = 1e6, overlap = 1)


###################################################
### code chunk number 7: BAFWin
###################################################
snp.b.win <- windowValues(x = sample.het.i$Bf,
                          positions = sample.het.i$position,
                          chromosomes = sample.het.i$chromosome,
                          window = 1e6, overlap = 1)


###################################################
### code chunk number 8: chrViewNoMut
###################################################

chromosome.view(baf.windows = snp.b.win[[1]], 
                ratio.windows = snp.r.win[[1]],
                 min.N.ratio = 1)


###################################################
### code chunk number 9: copynumber
###################################################
breaks <- find.breaks(sample.het.i, gamma = 20, kmin = 15, baf.thres = c(0, 0.5))
seg.i  <- segment.breaks(sample.i, breaks = breaks)


###################################################
### code chunk number 10: adjustSegs
###################################################

weights.snp   <- 150 + round((seg.i$end.pos - seg.i$start.pos)/1e6 , 0)
filter.size   <- (seg.i$end.pos - seg.i$start.pos) >= 10e6
avg.unlogR <- mean(sample.i$adjusted.ratio, na.rm = TRUE)


###################################################
### code chunk number 11: CPforSNP (eval = FALSE)
###################################################
## CPsnp.example <- baf.model.fit(Bf = seg.i$Bf[filter.size],
##                         depth.ratio = seg.i$depth.ratio[filter.size], 
##                         weight.ratio = weights.snp[filter.size],
##                         weight.Bf = weights.snp[filter.size],
##                         avg.depth.ratio = avg.unlogR,
##                         cellularity = seq(0.1,1,0.01), 
##                         ploidy = seq(1,7,0.1),
##                         priors.table = data.frame(CN = 2, value = 2))


###################################################
### code chunk number 12: loadCPforSNP
###################################################
data(CPsnp.example)


###################################################
### code chunk number 13: getCPparamSNP
###################################################
cint <- get.ci(CPsnp.example)

cellularity <- cint$max.cellularity
ploidy   <- cint$max.ploidy


###################################################
### code chunk number 14: CPplotSNP
###################################################
cp.plot(CPsnp.example)
cp.plot.contours(CPsnp.example, add = TRUE)


###################################################
### code chunk number 15: bafBayesSNP
###################################################
snp.seg.cn <- baf.bayes(Bf = seg.i$Bf,
                        depth.ratio = seg.i$depth.ratio, 
                        avg.depth.ratio = avg.unlogR,
                        cellularity = cellularity,
                        weight.ratio = 2 * 300,
                        weight.Bf = 300, ratio.priority = FALSE,
                        ploidy = ploidy, CNt.max = 10)

segmented.snp <- cbind(seg.i, snp.seg.cn)

head(segmented.snp[segmented.snp$chromosome == 1, ])


###################################################
### code chunk number 16: chrViewNoMutCP1
###################################################

chromosome.view(baf.windows = snp.b.win[[1]], 
                ratio.windows = snp.r.win[[1]],  min.N.ratio = 1,
                segments = segmented.snp[segmented.snp$chromosome == "1", ],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.unlogR, main = "1")


###################################################
### code chunk number 17: chrViewNoMutCP16
###################################################
chromosome.view(baf.windows = snp.b.win[[16]], 
                ratio.windows = snp.r.win[[16]],  min.N.ratio = 1,
                segments = segmented.snp[segmented.snp$chromosome == "16", ],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.unlogR, main = "16")


###################################################
### code chunk number 18: genomeViewCNtSNP
###################################################

genome.view(seg.cn = segmented.snp, info.type = "CNt")
legend("bottomright", bty="n", c("Tumor copy number"),col = c("red"), 
       inset = c(0, -0.4), pch=15, xpd = TRUE)


###################################################
### code chunk number 19: genomeViewABSNP
###################################################
genome.view(seg.cn = segmented.snp, info.type = "AB")
legend("bottomright", bty = "n", c("A-allele","B-allele"), col= c("red", "blue"), 
       inset = c(0, -0.45), pch = 15, xpd = TRUE)


