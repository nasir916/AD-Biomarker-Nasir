
## 5_sensitivity_analyses.R — sensitivity and secondary analyses (labelled as such in the manuscript)
##  A robustness grid (alpha, gene cap, filter size) and alternative scalings, each with full nested CV + validation
##  B frozen model applied to control-referenced and fixed-reference (single-sample) validation data
##  C leave-one-batch-out CV in discovery
##  D development-uncertainty bootstrap (500 pipeline refits applied to validation) + bootstrap CIs of validation metrics
##  E adjustment for age, sex and blood cell-type scores; threshold sensitivity
##  F specificity: cross-fitted scores for FTD-spectrum (GSE140830, duplicates of discovery arrays removed) and MCI
source("R/pipeline_functions.R"); library(parallel)
set.seed(20261025); dir.create("results/sensitivity", showWarnings = FALSE)
NCORE <- 16; R <- 10; K <- 5
sc <- c("cohortz","ctlref","fixedref")
fn <- function(c, s) if (s == "cohortz") sprintf("data/%s_z.csv.gz", c) else sprintf("data/%s_%s.csv.gz", c, s)  # cohortz = original analysis files
DAT <- lapply(setNames(sc, sc), function(s) list(
  disc = read_cohort(fn("GSE140829", s)), v61 = read_cohort(fn("GSE63061", s)), v60 = read_cohort(fn("GSE63060", s))))
for (s in sc) for (k in names(DAT[[s]])) stopifnot(identical(DAT[[s]][[k]]$id, DAT$cohortz[[k]]$id))
genes <- setdiff(names(DAT$cohortz$disc), META)
XM <- lapply(DAT, function(L) lapply(L, function(d) as.matrix(d[, genes])))
for (s in sc) for (k in names(DAT[[s]])) DAT[[s]][[k]] <- DAT[[s]][[k]][, intersect(META, names(DAT[[s]][[k]]))]
gc()
dup <- read.csv("results/GSE140830_duplicates_of_GSE140829.csv")
ftd <- read_cohort("data/GSE140830_fixedref.csv.gz"); ftd <- ftd[!(ftd$id %in% dup$id[dup$dup_of != ""]), ]
mci <- read_cohort("data/GSE140829_MCI_fixedref.csv.gz")
EXT <- rbind(data.frame(grp = ifelse(ftd$dx == "Control", "Independent control (GSE140830)", ftd$dx), as.data.frame(as.matrix(ftd[, genes])), check.names = FALSE),
             data.frame(grp = "MCI (GSE140829)", as.data.frame(as.matrix(mci[, genes])), check.names = FALSE))
XE <- as.matrix(EXT[, genes])
frozen <- readRDS("results/frozen_model.rds")

configs <- list(
  base = list(sc = "cohortz", alpha = 0.5, cap = 20, top = 200),
  alpha0.25 = list(sc = "cohortz", alpha = 0.25, cap = 20, top = 200),
  alpha0.75 = list(sc = "cohortz", alpha = 0.75, cap = 20, top = 200),
  lasso = list(sc = "cohortz", alpha = 1, cap = 20, top = 200),
  cap10 = list(sc = "cohortz", alpha = 0.5, cap = 10, top = 200),
  cap15 = list(sc = "cohortz", alpha = 0.5, cap = 15, top = 200),
  cap30 = list(sc = "cohortz", alpha = 0.5, cap = 30, top = 200),
  nocap = list(sc = "cohortz", alpha = 0.5, cap = Inf, top = 200),
  top100 = list(sc = "cohortz", alpha = 0.5, cap = 20, top = 100),
  top500 = list(sc = "cohortz", alpha = 0.5, cap = 20, top = 500),
  ctlref = list(sc = "ctlref", alpha = 0.5, cap = 20, top = 200),
  fixedref = list(sc = "fixedref", alpha = 0.5, cap = 20, top = 200))
tasks <- list()
for (cn in names(configs)) for (r in 1:R) { f <- strat_folds(DAT$cohortz$disc$AD, K, 1000 + r); for (k in 1:K) tasks[[length(tasks) + 1]] <- list(cn = cn, r = r, k = k, te = which(f == k)) }
for (cn in names(configs)) tasks[[length(tasks) + 1]] <- list(cn = cn, r = 0, k = 0, te = integer(0))   # final fit
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261025)
clusterExport(cl, c("DAT", "XM", "XE", "configs", "genes"))
clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
res <- parLapply(cl, tasks, function(t) {
  cf <- configs[[t$cn]]; d <- DAT[[cf$sc]]$disc; X <- XM[[cf$sc]]$disc
  tr <- setdiff(seq_len(nrow(d)), t$te)
  m <- pipeline(d, X, tr, cf$top, cf$alpha, cf$cap)
  out <- list(cn = t$cn, r = t$r, k = t$k, te = t$te, genes = m$genes, penalty = m$penalty)
  if (length(t$te)) out$lp_te <- lp_new(m, X[t$te, , drop = FALSE])
  if (t$cn == "fixedref" && length(t$te)) out$lp_ext <- lp_new(m, XE)
  if (t$r == 0) { out$lp61 <- lp_new(m, XM[[cf$sc]]$v61); out$lp60 <- lp_new(m, XM[[cf$sc]]$v60); out$model <- m }
  out })
## ---- A: summarise grid ----
y <- DAT$cohortz$disc$AD; y61 <- DAT$cohortz$v61$AD; y60 <- DAT$cohortz$v60$AD
G0 <- frozen$genes; grid <- list(); ext_scores <- matrix(0, nrow(EXT), 1); n_ext <- 0
for (cn in names(configs)) {
  rr <- Filter(function(x) x$cn == cn && x$r > 0, res); oof <- matrix(NA, length(y), R)
  for (x in rr) oof[x$te, x$r] <- x$lp_te
  cva <- sapply(1:R, function(r) as.numeric(auc(y, oof[, r], quiet = TRUE, direction = "<")))
  fin <- Filter(function(x) x$cn == cn && x$r == 0, res)[[1]]
  a61 <- auc_ci(y61, fin$lp61); a60 <- auc_ci(y60, fin$lp60); c61 <- cal(y61, fin$lp61); c60 <- cal(y60, fin$lp60)
  ng <- sapply(rr, function(x) length(x$genes))
  grid[[cn]] <- data.frame(config = cn, scaling = configs[[cn]]$sc, alpha = configs[[cn]]$alpha, cap = configs[[cn]]$cap, top = configs[[cn]]$top,
    cv_auc_mean = mean(cva), cv_auc_min = min(cva), cv_auc_max = max(cva), median_genes_fold = median(ng),
    final_genes = length(fin$genes), overlap_with_frozen17 = length(intersect(fin$genes, G0)),
    AUC61 = a61[2], AUC61_lo = a61[1], AUC61_hi = a61[3], slope61 = c61["slope"], int61 = c61["intercept"],
    AUC60 = a60[2], AUC60_lo = a60[1], AUC60_hi = a60[3], slope60 = c60["slope"], int60 = c60["intercept"])
  if (cn == "fixedref") {
    E <- do.call(cbind, lapply(rr, function(x) x$lp_ext)); EXT$lp_cf <- rowMeans(E)
    disc_cf <- data.frame(id = DAT$fixedref$disc$id, AD = y, lp_cf = rowMeans(oof)) }
  if (cn == "base") write.csv(data.frame(config = cn, gene = unlist(lapply(rr, `[[`, "genes"))), "results/sensitivity/base_rerun_selected_genes.csv", row.names = FALSE)
}
grid <- do.call(rbind, grid); rownames(grid) <- NULL
write.csv(grid, "results/sensitivity/A_robustness_grid.csv", row.names = FALSE); print(grid[, c(1, 6, 10, 11, 12, 13, 18)], digits = 3)
## ---- B: frozen model on alternative scalings of validation data ----
B <- list()
for (s in sc) for (v in c("v61", "v60")) {
  yy <- DAT[[s]][[v]]$AD; lp <- lp_new(frozen$final, XM[[s]][[v]]); a <- auc_ci(yy, lp); cc <- cal(yy, lp)
  pcls <- as.integer(plogis(lp) >= frozen$threshold)
  B[[paste(s, v)]] <- data.frame(scaling = s, cohort = ifelse(v == "v61", "GSE63061", "GSE63060"), AUC = a[2], lo = a[1], hi = a[3],
    cal_intercept = cc["intercept"], cal_slope = cc["slope"], sens = mean(pcls[yy == 1]), spec = mean(1 - pcls[yy == 0])) }
B <- do.call(rbind, B); rownames(B) <- NULL; write.csv(B, "results/sensitivity/B_frozen_model_scalings.csv", row.names = FALSE); print(B, digits = 3)
## ---- C: leave-one-batch-out ----
d <- DAT$cohortz$disc; X <- XM$cohortz$disc; bl <- levels(d$batch)
clusterExport(cl, c("d", "X"))
lobo <- parLapply(cl, bl, function(b) { tr <- which(d$batch != b); te <- which(d$batch == b)
  m <- pipeline(d, X, tr); list(b = b, te = te, lp = lp_new(m, X[te, , drop = FALSE]), genes = m$genes) })
lp_lobo <- rep(NA, nrow(d)); for (x in lobo) lp_lobo[x$te] <- x$lp
pb <- do.call(rbind, lapply(lobo, function(x) { yy <- y[x$te]; data.frame(batch = x$b, n = length(yy), n_AD = sum(yy), n_CN = sum(1 - yy),
  AUC = if (sum(yy) >= 5 && sum(1 - yy) >= 5) as.numeric(auc(yy, x$lp, quiet = TRUE, direction = "<")) else NA,
  overlap_with_frozen17 = length(intersect(x$genes, G0))) }))
pooled <- auc_ci(y, lp_lobo)
## batch-stratified AUC: probability a random AD case scores above a random CN from the same batch
num <- 0; den <- 0
for (b in bl) { i1 <- which(d$batch == b & y == 1); i0 <- which(d$batch == b & y == 0)
  if (length(i1) && length(i0)) { o <- outer(lp_lobo[i1], lp_lobo[i0], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) } }
write.csv(pb, "results/sensitivity/C_lobo_per_batch.csv", row.names = FALSE)
lobo_sum <- data.frame(metric = c("LOBO pooled AUC", "LOBO pooled AUC lo", "LOBO pooled AUC hi", "LOBO within-batch (stratified) AUC", "Batches", "Batches with >=5 AD and >=5 CN", "Median per-batch AUC (eligible)"),
  value = c(pooled[2], pooled[1], pooled[3], num / den, length(bl), sum(!is.na(pb$AUC)), median(pb$AUC, na.rm = TRUE)))
write.csv(lobo_sum, "results/sensitivity/C_lobo_summary.csv", row.names = FALSE); print(lobo_sum, digits = 3)

## ---- C2: LOBO with within-batch standardisation (mirrors within-cohort scaling used for external cohorts) ----
Xb <- X
for (b in levels(d$batch)) { i <- which(d$batch == b); if (length(i) > 2) Xb[i, ] <- scale(X[i, ]) else Xb[i, ] <- sweep(X[i, , drop = FALSE], 2, colMeans(X[i, , drop = FALSE])) }
Xb[is.na(Xb)] <- 0
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261025)
clusterExport(cl, c("d", "Xb")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
lobo2 <- parLapply(cl, levels(d$batch), function(b) { tr <- which(d$batch != b); te <- which(d$batch == b)
  m <- pipeline(d, Xb, tr); list(b = b, te = te, lp = lp_new(m, Xb[te, , drop = FALSE]), genes = m$genes) })
stopCluster(cl)
lp2 <- rep(NA, nrow(d)); for (x in lobo2) lp2[x$te] <- x$lp
num <- 0; den <- 0
for (b in levels(d$batch)) { i1 <- which(d$batch == b & d$AD == 1); i0 <- which(d$batch == b & d$AD == 0)
  if (length(i1) && length(i0)) { o <- outer(lp2[i1], lp2[i0], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) } }
pc <- auc_ci(d$AD, lp2)
v2 <- data.frame(metric = c("LOBO (within-batch standardised) pooled AUC", "lo", "hi", "LOBO (within-batch standardised) within-batch AUC", "Median overlap with frozen 17 genes"),
                 value = c(pc[2], pc[1], pc[3], num / den, median(sapply(lobo2, function(x) length(intersect(x$genes, G0))))))
write.csv(v2, "results/sensitivity/C_lobo_batchstandardised_summary.csv", row.names = FALSE)

## ---- D: development-uncertainty bootstrap ----
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261025)
clusterExport(cl, c("DAT", "XM", "d", "X")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
NB <- 500
dev <- parLapply(cl, 1:NB, function(b) { set.seed(5000 + b); bi <- sample(seq_len(nrow(d)), replace = TRUE)
  m <- tryCatch(pipeline(d, X, bi), error = function(e) NULL); if (is.null(m)) return(c(NA, NA, NA))
  c(as.numeric(auc(DAT$cohortz$v61$AD, lp_new(m, XM$cohortz$v61), quiet = TRUE, direction = "<")),
    as.numeric(auc(DAT$cohortz$v60$AD, lp_new(m, XM$cohortz$v60), quiet = TRUE, direction = "<")), length(m$genes)) })
stopCluster(cl)
dev <- do.call(rbind, dev); colnames(dev) <- c("AUC61", "AUC60", "n_genes"); write.csv(dev, "results/sensitivity/D_development_bootstrap_draws.csv", row.names = FALSE)
## bootstrap CIs of frozen-model validation metrics (resampling validation participants)
fv <- read.csv("results/validation_predictions.csv")
bootv <- list()
for (cn in c("GSE63061", "GSE63060")) {
  v <- fv[fv$cohort == cn, ]; vv <- DAT$cohortz[[ifelse(cn == "GSE63061", "v61", "v60")]]; v$age <- vv$age[match(v$id, vv$id)]; v$sex <- vv$sex[match(v$id, vv$id)]
  stat <- function(ix) { w <- v[ix, ]; c(as.numeric(auc(w$AD, w$lp, quiet = TRUE, direction = "<")), NA, NA,   # calibration not bootstrapped here: outcome-stratified resampling fixes prevalence
    { a0 <- glm(AD ~ age + sex, binomial, w); a1 <- glm(AD ~ age + sex + lp, binomial, w)
      as.numeric(auc(w$AD, fitted(a1), quiet = TRUE, direction = "<")) - as.numeric(auc(w$AD, fitted(a0), quiet = TRUE, direction = "<")) }) }
  set.seed(7); S <- t(replicate(2000, { ix <- c(sample(which(v$AD == 1), replace = TRUE), sample(which(v$AD == 0), replace = TRUE)); stat(ix) }))
  est <- stat(seq_len(nrow(v))); tp <- sum(v$p >= frozen$threshold & v$AD == 1); tn <- sum(v$p < frozen$threshold & v$AD == 0)
  ws <- wilson(tp, sum(v$AD)); wp <- wilson(tn, sum(1 - v$AD))
  bootv[[cn]] <- data.frame(cohort = cn, metric = c("AUC", "Calibration intercept", "Calibration slope", "Delta AUC (age+sex -> +signature, in-cohort)", "Sensitivity", "Specificity"),
    estimate = c(est, tp / sum(v$AD), tn / sum(1 - v$AD)),
    lo = c(apply(S, 2, quantile, 0.025, na.rm = TRUE), ws[1], wp[1]), hi = c(apply(S, 2, quantile, 0.975, na.rm = TRUE), ws[2], wp[2]),
    method = c(rep("bootstrap 2000 (stratified)", 4), rep("Wilson", 2)),
    dev_lo = c(quantile(dev[, ifelse(cn == "GSE63061", "AUC61", "AUC60")], 0.025, na.rm = TRUE), rep(NA, 5)),
    dev_hi = c(quantile(dev[, ifelse(cn == "GSE63061", "AUC61", "AUC60")], 0.975, na.rm = TRUE), rep(NA, 5)),
    dev_median = c(median(dev[, ifelse(cn == "GSE63061", "AUC61", "AUC60")], na.rm = TRUE), rep(NA, 5))) }
bootv <- do.call(rbind, bootv); write.csv(bootv, "results/sensitivity/D_validation_uncertainty.csv", row.names = FALSE); print(bootv, digits = 3)
## ---- E: adjustment for age, sex, cell-type scores; threshold sensitivity ----
mk <- list(Neutrophil = c("FCGR3B","CSF3R","CXCR2","FPR1"), Monocyte = c("CD14","LYZ","CSF1R","FCN1"), Bcell = c("MS4A1","CD19","CD79A","BANK1"),
           Tcell = c("CD3E","CD3D","CD2","LCK"), NK = c("GNLY","NKG7","KLRD1","PRF1"))
oofd <- read.csv("results/discovery_oof_predictions.csv")
adj <- list()
for (cn in c("GSE140829", "GSE63061", "GSE63060")) {
  key <- switch(cn, GSE140829 = "disc", GSE63061 = "v61", GSE63060 = "v60"); w <- DAT$cohortz[[key]]; Xk <- XM$cohortz[[key]]
  for (k in names(mk)) w[[k]] <- rowMeans(Xk[, intersect(mk[[k]], colnames(Xk)), drop = FALSE])
  w$lp <- if (cn == "GSE140829") oofd$oof_lp_mean[match(w$id, oofd$id)] else fv$lp[fv$cohort == cn][match(w$id, fv$id[fv$cohort == cn])]
  w$lpz <- as.numeric(scale(w$lp))
  base_f <- AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK
  if (cn == "GSE140829") base_f <- update(base_f, . ~ . + batch)
  m0 <- glm(base_f, binomial, w); m1 <- update(m0, . ~ . + lpz); m_u <- glm(AD ~ lpz, binomial, w)
  lr <- anova(m0, m1, test = "LRT"); b <- coef(summary(m1))["lpz", ]
  r0 <- roc(w$AD, fitted(m0), quiet = TRUE, direction = "<"); r1 <- roc(w$AD, fitted(m1), quiet = TRUE, direction = "<")
  adj[[cn]] <- data.frame(cohort = cn, lp_source = ifelse(cn == "GSE140829", "out-of-fold", "frozen"),
    OR_unadjusted = exp(coef(m_u)["lpz"]), OR_adjusted = exp(b[1]), OR_adj_lo = exp(b[1] - 1.96 * b[2]), OR_adj_hi = exp(b[1] + 1.96 * b[2]),
    LR_chisq = lr$Deviance[2], LR_p = lr$`Pr(>Chi)`[2], AUC_covariates = as.numeric(auc(r0)), AUC_cov_plus_sig = as.numeric(auc(r1)),
    dAUC_p = roc.test(r1, r0, paired = TRUE)$p.value, adjusted_for = paste(all.vars(base_f)[-1], collapse = "+")) }
adj <- do.call(rbind, adj); rownames(adj) <- NULL; write.csv(adj, "results/sensitivity/E_adjusted_cell_composition.csv", row.names = FALSE); print(adj, digits = 3)
thr <- list()
for (cn in c("GSE63061", "GSE63060")) { v <- fv[fv$cohort == cn, ]
  for (t in c(0.35, 0.40, 0.451, 0.50, 0.55, 0.60)) thr[[paste(cn, t)]] <- data.frame(cohort = cn, threshold = t, frozen = t == 0.451,
    sens = mean(v$p[v$AD == 1] >= t), spec = mean(v$p[v$AD == 0] < t), PPV = mean(v$AD[v$p >= t]), NPV = mean(1 - v$AD[v$p < t])) }
thr <- do.call(rbind, thr); rownames(thr) <- NULL; write.csv(thr, "results/sensitivity/E_threshold_sensitivity.csv", row.names = FALSE)
## ---- F: specificity (cross-fitted fold models, fixed-reference scaling) ----
allsc <- rbind(data.frame(grp = ifelse(disc_cf$AD == 1, "AD (GSE140829)", "Control (GSE140829)"), lp = disc_cf$lp_cf),
               data.frame(grp = EXT$grp, lp = EXT$lp_cf))
write.csv(allsc, "results/sensitivity/F_crossfitted_scores_by_group.csv", row.names = FALSE)
cmp <- function(a, b) { s1 <- allsc$lp[allsc$grp == a]; s0 <- allsc$lp[allsc$grp == b]
  ci <- auc_ci(c(rep(1, length(s1)), rep(0, length(s0))), c(s1, s0)); data.frame(case = a, reference = b, n_case = length(s1), n_ref = length(s0), AUC = ci[2], lo = ci[1], hi = ci[3]) }
grps <- c("bvFTD", "PSP", "nfvPPA", "svPPA", "CBS", "MCI (GSE140829)")
spec <- rbind(cmp("AD (GSE140829)", "Control (GSE140829)"),
  do.call(rbind, lapply(grps, function(g) cmp(g, "Control (GSE140829)"))),
  do.call(rbind, lapply(grps[1:5], function(g) cmp("AD (GSE140829)", g))),
  cmp("Independent control (GSE140830)", "Control (GSE140829)"))
## (group AUCs are reported from the matched-ensemble construction in 7_secondary_analyses.R)
writeLines(capture.output(sessionInfo()), "results/sensitivity/sessionInfo_R.txt")
