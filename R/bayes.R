mufreq.dbinom <- function(mufreq, mufreq.model, depth.t,
    seq.errors = 0.01, ...) {
    mufreq.model[mufreq.model == 0] <- seq.errors
    n.success <- round(mufreq * depth.t, 0)
    dbinom(x = n.success, size = depth.t, prob = mufreq.model, ...)
}

mufreq.dpois <- function(mufreq, mufreq.model, depth.t,
    seq.errors = 0.01, ...) {
    mufreq.model[mufreq.model == 0] <- seq.errors
    n.success <- round(mufreq * depth.t, 0)
    dpois( x = n.success, lambda = mufreq.model * depth.t, ...)
}

baf.dbinom <- function(baf, baf.model, depth.t, ...) {
    n.success <- round(baf * depth.t, 0)
    dbinom( x = n.success, size = depth.t, prob = baf.model, ...)
}

baf.dpois <- function(baf, baf.model, depth.t, ...) {
    n.success <- round(baf * depth.t, 0)
    dpois( x = n.success, lambda = baf.model * depth.t, ...)
}

depth.ratio.dbinom <- function(size, depth.ratio, depth.ratio.model, ...) {
    n.success <- round(size * (depth.ratio / (1 + depth.ratio)), 0)
    prob <- depth.ratio.model / (1 + depth.ratio.model)
    dbinom( x = n.success, size = size, prob = prob, ...)
}

depth.ratio.dpois <- function(size, depth.ratio, depth.ratio.model, ...) {
    x <- round(size * depth.ratio, 0)
    dpois( x = x, lambda = depth.ratio.model * size, ...)
}



baf.bayes <- function(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
    sd.Bf = 0.1, sd.ratio = 0.5, weight.Bf = 1, weight.ratio = 1, CNt.min = 0,
    CNt.max = 7, CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max,
    value = 1), ratio.priority = FALSE) {

    baf.tab <- data.frame(Bf = Bf, ratio = log(depth.ratio),
        sd.Bf = sd.Bf, sd.ratio = sd.ratio, weight.Bf = weight.Bf,
        weight.ratio = weight.ratio)
    baf_types <- baf.types.matrix(CNt.min = CNt.min,
        CNt.max = CNt.max, CNn = CNn)
    model_baf <- baf.model.points(cellularity = cellularity,
        ploidy = ploidy, baf_types = baf_types,
        avg.depth.ratio = avg.depth.ratio)
    model.pts <- cbind(CNt =  baf_types$CNt, A = baf_types$CNt - baf_types$B,
        B = baf_types$B, model_baf)

    rows.x <- 1:nrow(baf.tab)
    priors <- rep(1, nrow(model.pts))
    for (i in 1:nrow(priors.table)) {
        priors[model.pts$CNt == priors.table$CN[i]] <- priors.table$value[i]
    }
    priors <- priors / sum(priors)

    bayes.fit <- function (x, mat, model.pts, priors, ratio.priority) {
        test.ratio <- log(model.pts$depth.ratio)
        test.baf <- model.pts$BAF
        min.offset <- 1e-323
        score.r <- dt2(sd = mat[x, ]$sd.ratio / sqrt(mat[x, ]$weight.ratio),
            mean = mat[x, ]$ratio, x = test.ratio, df = 5, log = TRUE)
        score.r <- score.r + log(priors)
        if (!is.na(mat[x, ]$Bf) & !is.na(mat[x, ]$sd.Bf /
            sqrt(mat[x, ]$weight.Bf))) {
            score.b <- dt2(mean = mat[x, ]$Bf, sd = mat[x, ]$sd.Bf /
                sqrt(mat[x, ]$weight.Bf), x = test.baf, df = 5, log = TRUE)
            post.model <- score.r + score.b
        } else {
            post.model <- score.r
        }
        if (ratio.priority == FALSE) {
            max.lik <-  which.max(post.model)
            max.post <- c(as.numeric(model.pts[max.lik, 1:3]),
                post.model[max.lik])
        } else {
            res.cn <- model.pts$CNt[which.max(score.r)]
            idx.pts <- model.pts$CNt == res.cn
            model.lik <- cbind(model.pts[idx.pts, 1:3], post.model[idx.pts])
            if (is.null(dim(model.lik))) {
                max.post <- model.lik
            } else {
                max.post   <- model.lik[which.max(model.lik[, 4]), ]
            }
        }
        if (is.na(mat[x, ]$Bf)) {
             max.post[2:3] <- NA
        }
        max.post
    }
    bafs.L <- mapply(FUN = bayes.fit, rows.x, MoreArgs = list(mat = baf.tab,
        model.pts = model.pts, priors = priors,
        ratio.priority = ratio.priority), SIMPLIFY = FALSE)
    bafs.L <- do.call(rbind, bafs.L)
    colnames(bafs.L) <- c("CNt", "A", "B", "LPP")
    bafs.L
}

baf.model.fit <- function(cellularity = seq(0.3, 1, by = 0.01),
    ploidy = seq(1, 7, by = 0.1), mc.cores = getOption("mc.cores", 2L), ...) {

    result <- expand.grid(ploidy = ploidy, cellularity = cellularity,
        KEEP.OUT.ATTRS = FALSE)

    fit.cp <- function(ii) {
        L.model <- baf.bayes(cellularity = result$cellularity[ii],
            ploidy = result$ploidy[ii], ...)
        sum(L.model[,4])
    }
    bayes.res <- pblapply(X = 1:nrow(result), FUN = fit.cp,
        cl = mc.cores)
    result$LPP <- unlist(bayes.res)
    z <- tapply(result$LPP, list(result$ploidy, result$cellularity), mean)
    x <- as.numeric(rownames(z))
    y <- as.numeric(colnames(z))
    max.lik <- max(result$LPP, na.rm = TRUE)
    LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
    znorm <- exp(z - LogSumLik)
    list(ploidy = x, cellularity = y, lpp = znorm)
}

mufreq.bayes <- function(mufreq, depth.ratio, cellularity,
    ploidy, avg.depth.ratio, weight.mufreq = 100, weight.ratio = 100,
    CNt.min = 1, CNt.max = 7, CNn = 2,
    priors.table = data.frame(CN = CNt.min:CNt.max, value = 1)) {
    mufreq.tab <- data.frame(F = mufreq, ratio = depth.ratio,
        weight.mufreq = weight.mufreq, weight.ratio = weight.ratio)
    mufreq_types <- mufreq.types.matrix(CNt.min = CNt.min,
        CNt.max = CNt.max, CNn = CNn)
    model.pts <- mufreq.model.points(cellularity = cellularity,
        ploidy = ploidy, mufreq_types = mufreq_types,
        avg.depth.ratio = avg.depth.ratio)
    model.pts <- cbind(mufreq_types, model.pts)
    rows.x <- 1:nrow(mufreq.tab)
    priors <- rep(1, nrow(model.pts))
    for (i in 1:nrow(priors.table)) {
        priors[model.pts$CNt ==
            priors.table$CN[i]] <- priors.table$value[i]
    }
    priors <- priors / sum(priors)
    bayes.fit <- function (x, mat, model.pts, priors) {
        test.ratio <- model.pts$depth.ratio
        test.mufrq <- model.pts$mufreqs
        min.offset <- 1e-323
        score.r <- depth.ratio.dbinom(size = mat[x, ]$weight.ratio,
            depth.ratio = mat[x, ]$ratio, test.ratio)
        score.m <- mufreq.dbinom(mufreq = mat[x, ]$F,
            depth.t = mat[x, ]$weight.mufreq, test.mufrq)
        score.r <- score.r * priors
        score.m <- score.m
        post.model <- score.r * score.m
        post.model[post.model == 0] <- min.offset
        res.cn <- model.pts$CNt[which.max(score.r)]
        idx.pts <- model.pts$CNt == res.cn
        model.lik <- cbind(model.pts[idx.pts, 1:3], log(post.model[idx.pts]))
        if (is.null(dim(model.lik))) {
            max.post <- model.lik
        } else {
            max.post <- model.lik[which.max(model.lik[,4]),]
        }
        max.post
    }
    types.L <- mapply(FUN = bayes.fit, rows.x, MoreArgs = list(
        mat = mufreq.tab, model.pts = model.pts,
        priors = priors), SIMPLIFY = FALSE)
   types.L <- do.call(rbind, types.L)
   colnames(types.L) <- c("CNn", "CNt", "Mt", "LPP")
   types.L
}

mufreq.model.fit <- function(cellularity = seq(0.3, 1, by = 0.01),
    ploidy = seq(1, 7, by = 0.1),  mc.cores = getOption("mc.cores", 2L), ...) {
   result <- expand.grid(ploidy = ploidy, cellularity = cellularity,
        KEEP.OUT.ATTRS = FALSE)
    fit.cp <- function(ii) {
        L.model <- mufreq.bayes(cellularity = result$cellularity[ii],
            ploidy = result$ploidy[ii], ...)
        sum(L.model[, 4])
    }
    bayes.res <- pblapply(X = 1:nrow(result), FUN = fit.cp,
        cl = mc.cores)
    result$LPP <- unlist(bayes.res)
    z <- tapply(result$LPP, list(result$ploidy, result$cellularity), mean)
    x <- as.numeric(rownames(z))
    y <- as.numeric(colnames(z))
    max.lik <- max(result$LPP, na.rm = TRUE)
    LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
    znorm <- exp(z - LogSumLik)
    list(ploidy = x, cellularity = y, lpp = znorm)
}
