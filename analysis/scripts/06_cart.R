#!/usr/bin/env Rscript
# Step 6: CART analyses (rpart = Breiman CART, same algorithm as sklearn trees).
#
#  mode "faithful": CART-CELLS settings (Carluer et al. 2024, github CarluerJB/
#     CART_CELLS, PARAM_CART_DEFAULT): classification tree on raw counts,
#     features = GR-TF counts, max_depth 3, min_samples_split = min_samples_leaf
#     = 10% of cells, Gini, no pruning; 80/20 split (seed 42) for accuracy;
#     Mann-Whitney on the TF1 split; TF1 x TF2 and TF1 x TF3 two-way ANOVA.
#     Run independently in Exp1 and Exp2 -> reproducibility of TF1/TF2 calls.
#  mode "grid": regression CART on covariate-residualised log expression with
#     binary GR-TF presence (>=1 read) as features; grid over maxdepth,
#     minbucket and cp; each tree trained in one experiment predicts the other
#     (Exp1 -> Exp2 and Exp2 -> Exp1). Score = per-gene correlation between
#     prediction (partialled on test covariates) and observed expression.
#  mode "tuned": the best grid setting on all genes, both directions, plus the
#     CART-CELLS setting for comparison; saves per-gene transfer r and trees.
suppressPackageStartupMessages({library(rpart); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
mode <- ifelse(length(args) >= 1, args[1], "faithful")
ncores <- as.integer(Sys.getenv("NCORES", "10"))
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
fr <- file.path(base, "data", "forR"); out <- file.path(base, "results", "cart"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
genes <- readLines(file.path(fr, "genes.txt")); tfs <- readLines(file.path(fr, "tfs.txt"))
G <- length(genes); K <- length(tfs)
ncell <- function(e) file.size(file.path(fr, paste0(e, "_TF.i32"))) / 4 / K
rd <- function(e, what, type, ncol) {
  n <- ncell(e)
  matrix(readBin(file.path(fr, paste0(e, "_", what)), what = type, n = n * ncol, size = 4), n, ncol)
}
load_exp <- function(e, raw = FALSE) {
  n <- ncell(e)
  TF <- rd(e, "TF.i32", "integer", K); colnames(TF) <- make.names(tfs)
  nz <- file.size(file.path(fr, paste0(e, "_Z.f32"))) / 4 / n
  Z <- rd(e, "Z.f32", "numeric", nz)
  lst <- list(n = n, TF = TF, Z = Z, Yres = rd(e, "Yres.f32", "numeric", G))
  if (raw) lst$Yraw <- rd(e, "Yraw.i32", "integer", G)
  lst
}

# ---------------------------------------------------------------- faithful
cart_cells_gene <- function(y, TF, n) {
  ms <- ceiling(0.10 * n)
  # cp = -1: grow like sklearn (split whenever Gini decreases; rpart cp counts misclassification)
  ctrl <- rpart.control(minsplit = ms, minbucket = ms, maxdepth = 3, cp = -1, xval = 0,
                        maxcompete = 0, maxsurrogate = 0, usesurrogate = 0)
  df <- data.frame(y = factor(y), TF)
  fit <- tryCatch(rpart(y ~ ., data = df, method = "class", parms = list(split = "gini"), control = ctrl),
                  error = function(e) NULL)
  res <- list(TF1 = NA, TF2 = NA, TF3 = NA, lim1 = NA, lim2 = NA, lim3 = NA, p1 = NA, pN2 = NA, pN3 = NA,
              acc = NA, zero = mean(y == 0), mean1p = NA, mean1m = NA)
  if (is.null(fit) || nrow(fit$frame) < 3) return(res)
  fr_ <- fit$frame; nodes <- as.integer(rownames(fr_))
  sp <- fit$splits; spl_nodes <- nodes[fr_$var != "<leaf>"]
  getsplit <- function(nd) {
    i <- which(spl_nodes == nd); if (!length(i)) return(NULL)
    list(var = as.character(fr_$var[nodes == nd]), thr = sp[i, "index"])
  }
  s1 <- getsplit(1L); s2 <- getsplit(2L); s3 <- getsplit(3L)
  # CART-CELLS: TF2 = split under "TF1 <= lim" ; TF3 = split under "TF1 > lim".
  # rpart left child (node 2) is "< thr" when ncat = -1, ">= thr" when ncat = +1.
  left_is_low <- sp[which(spl_nodes == 1L), "ncat"] < 0
  sLow <- if (left_is_low) s2 else s3; sHigh <- if (left_is_low) s3 else s2
  res$TF1 <- s1$var; res$lim1 <- s1$thr
  hi1 <- TF[, s1$var] >= s1$thr
  ly <- log10(y); ly[!is.finite(ly)] <- 0
  res$p1 <- tryCatch(wilcox.test(ly[hi1], ly[!hi1])$p.value, error = function(e) NA)
  res$mean1p <- mean(y[hi1]); res$mean1m <- mean(y[!hi1])
  two_way <- function(s) {
    if (is.null(s)) return(NA)
    f1 <- factor(hi1); f2 <- factor(TF[, s$var] >= s$thr)
    if (nlevels(f2) < 2) return(NA)
    a <- anova(lm(y ~ f1 * f2)); p <- a[["Pr(>F)"]][3]; ifelse(is.null(p), NA, p)
  }
  if (!is.null(sLow)) { res$TF2 <- sLow$var; res$lim2 <- sLow$thr; res$pN2 <- two_way(sLow) }
  if (!is.null(sHigh)) { res$TF3 <- sHigh$var; res$lim3 <- sHigh$thr; res$pN3 <- two_way(sHigh) }
  set.seed(42); tr <- sample.int(n, floor(0.8 * n))
  fe <- tryCatch(rpart(y ~ ., data = df[tr, ], method = "class", parms = list(split = "gini"), control = ctrl),
                 error = function(e) NULL)
  if (!is.null(fe)) res$acc <- mean(as.character(predict(fe, df[-tr, ], type = "class")) == as.character(df$y[-tr]))
  res
}

if (mode == "faithful") {
  for (e in c("Exp1", "Exp2")) {
    d <- load_exp(e, raw = TRUE)
    t0 <- Sys.time()
    rs <- mclapply(seq_len(G), function(g) cart_cells_gene(d$Yraw[, g], d$TF, d$n), mc.cores = ncores)
    tab <- do.call(rbind, lapply(rs, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
    tab <- cbind(gene = genes, tab)
    write.table(tab, gzfile(file.path(out, paste0("cart_cells_faithful_", e, ".tsv.gz"))), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(e, "done in", format(Sys.time() - t0), "\n")
  }
}

# ---------------------------------------------------------------- transfer
partial_cor <- function(P, Yt, Zt) {
  # correlation of predictions with observed residualised expression, after
  # partialling the test covariates out of the predictions
  q <- qr(cbind(1, Zt)); Pr <- qr.resid(q, P)
  num <- colSums(Pr * Yt); den <- sqrt(colSums(Pr^2) * colSums(Yt^2))
  r <- num / den; r[!is.finite(r)] <- 0; r
}
fit_predict <- function(gidx, tr, te, depth, mb, cp, ytrain = "Yres") {
  Xtr <- as.data.frame((tr$TF >= 1) * 1); Xte <- as.data.frame((te$TF >= 1) * 1)
  ctrl <- rpart.control(maxdepth = depth, minbucket = mb, minsplit = 2 * mb, cp = cp, xval = 0,
                        maxcompete = 0, maxsurrogate = 0, usesurrogate = 0)
  P <- do.call(cbind, mclapply(gidx, function(g) {
    f <- rpart(y ~ ., data = cbind(y = tr[[ytrain]][, g], Xtr), method = "anova", control = ctrl)
    as.numeric(predict(f, Xte))
  }, mc.cores = ncores))
  partial_cor(P, te$Yres[, gidx, drop = FALSE], te$Z)
}

if (mode == "grid") {
  d1 <- load_exp("Exp1"); d2 <- load_exp("Exp2")
  set.seed(1); gidx <- sort(sample.int(G, as.integer(Sys.getenv("NGRID", "1500"))))
  grid <- expand.grid(depth = c(1, 2, 3, 4, 5), mb = c(25, 50, 100, 200, -0.10), cp = c(0, 0.001, 0.005, 0.01))
  res <- list()
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    for (dir in c("1to2", "2to1")) {
      tr <- if (dir == "1to2") d1 else d2; te <- if (dir == "1to2") d2 else d1
      mb <- if (g$mb < 0) ceiling(-g$mb * tr$n) else g$mb
      r <- fit_predict(gidx, tr, te, g$depth, mb, g$cp)
      n <- te$n; z <- atanh(r) * sqrt(n - 3); p <- pnorm(z, lower.tail = FALSE)
      res[[length(res) + 1]] <- data.frame(depth = g$depth, minbucket = g$mb, cp = g$cp, direction = dir,
        mean_r = mean(r), median_r = median(r), n_sig_bh = sum(p.adjust(p, "BH") < 0.05), n_genes = length(r))
    }
    cat(sprintf("%d/%d depth=%d mb=%s cp=%.3f  r(1to2)=%.4f r(2to1)=%.4f sig=%d/%d\n", i, nrow(grid), g$depth,
                g$mb, g$cp, res[[length(res) - 1]]$mean_r, res[[length(res)]]$mean_r,
                res[[length(res) - 1]]$n_sig_bh, res[[length(res)]]$n_sig_bh))
    write.table(do.call(rbind, res), file.path(out, "cart_grid_transfer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

if (mode == "tuned") {
  d1 <- load_exp("Exp1"); d2 <- load_exp("Exp2")
  d1$Ylog <- rd("Exp1", "Ylog.f32", "numeric", G); d2$Ylog <- rd("Exp2", "Ylog.f32", "numeric", G)
  depth <- as.integer(Sys.getenv("DEPTH")); mb <- as.numeric(Sys.getenv("MB")); cp <- as.numeric(Sys.getenv("CP"))
  # name = c(depth, minbucket (<0: fraction of cells), cp, train on adjusted (1) or raw log (0))
  settings <- list(tuned = c(depth, mb, cp, 1), stump = c(1, mb, cp, 1), cart_cells_unadjusted = c(3, -0.10, 0, 0))
  allr <- data.frame(gene = genes)
  for (s in names(settings)) for (dir in c("1to2", "2to1")) {
    tr <- if (dir == "1to2") d1 else d2; te <- if (dir == "1to2") d2 else d1
    v <- settings[[s]]; m <- if (v[2] < 0) ceiling(-v[2] * tr$n) else v[2]
    allr[[paste0("r_", s, "_", dir)]] <- fit_predict(seq_len(G), tr, te, v[1], m, v[3], ifelse(v[4] == 1, "Yres", "Ylog"))
    cat(s, dir, "mean r", mean(allr[[paste0("r_", s, "_", dir)]]), "\n")
  }
  mbn <- function(d) if (mb < 0) ceiling(-mb * d$n) else mb
  # tuned trees fitted within each experiment: root split + all split variables
  for (e in c("Exp1", "Exp2")) {
    d <- if (e == "Exp1") d1 else d2; X <- as.data.frame((d$TF >= 1) * 1)
    ctrl <- rpart.control(maxdepth = depth, minbucket = mbn(d), minsplit = 2 * mbn(d), cp = cp, xval = 0,
                          maxcompete = 0, maxsurrogate = 0, usesurrogate = 0)
    info <- mclapply(seq_len(G), function(g) {
      f <- rpart(y ~ ., data = cbind(y = d$Yres[, g], X), method = "anova", control = ctrl)
      v <- as.character(f$frame$var); v <- v[v != "<leaf>"]
      imp <- if (length(f$variable.importance)) f$variable.importance else numeric(0)
      c(root = ifelse(length(v), v[1], NA), vars = paste(unique(v), collapse = ","),
        n_splits = length(v), r2 = 1 - sum(residuals(f)^2) / sum((d$Yres[, g] - mean(d$Yres[, g]))^2),
        top_importance = ifelse(length(imp), names(imp)[1], NA))
    }, mc.cores = ncores)
    info <- do.call(rbind, info); colnames(info) <- paste0(colnames(info), "_", e)
    allr <- cbind(allr, info)
  }
  write.table(allr, gzfile(file.path(out, "cart_tuned_transfer_per_gene.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
}
