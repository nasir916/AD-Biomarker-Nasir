
## 10_additional_sensitivity_analyses.R — further sensitivity and secondary analyses (run after 5_ and 7_ in the same R session;
## run_all.py does this). Every analysis uses the discovery-developed model applied unchanged to the validation cohorts
## unless stated otherwise. Output: results/sensitivity/Q_*.csv ... Z_*.csv
if (!exists("bench2")) { source("5_sensitivity_analyses.R"); source("7_secondary_analyses.R") }
suppressPackageStartupMessages(library(parallel)); OUT <- "results/sensitivity/"
fm <- readRDS("results/frozen_model.rds"); G0 <- fm$genes
d <- DAT$cohortz$disc; X <- XM$cohortz$disc; y <- d$AD
V <- list(GSE63061 = list(d = DAT$cohortz$v61, X = XM$cohortz$v61), GSE63060 = list(d = DAT$cohortz$v60, X = XM$cohortz$v60))
auc1 <- function(yy, s) as.numeric(auc(yy, s, quiet = TRUE, direction = "<"))
pairauc <- function(s1, s0) { o <- outer(s1, s0, "-"); mean(o > 0) + 0.5 * mean(o == 0) }
mk <- list(Neutrophil = c("FCGR3B","CSF3R","CXCR2","FPR1"), Monocyte = c("CD14","LYZ","CSF1R","FCN1"), Bcell = c("MS4A1","CD19","CD79A","BANK1"),
           Tcell = c("CD3E","CD3D","CD2","LCK"), NK = c("GNLY","NKG7","KLRD1","PRF1"))
cscore <- function(Xk) sapply(mk, function(g) rowMeans(Xk[, intersect(g, colnames(Xk)), drop = FALSE]))

## ---- Q. fold-safe batch-residualised development (batch means estimated in training data only) ----
bres <- function(X, b, tr) {            # subtract training-set batch means; batches absent from training are centred on their own (label-free) mean
  Xr <- X; b <- as.character(b)
  for (bb in unique(b)) { i <- which(b == bb); itr <- intersect(i, tr)
    mu <- if (length(itr) >= 2) colMeans(X[itr, , drop = FALSE]) else colMeans(X[i, , drop = FALSE]); Xr[i, ] <- sweep(X[i, , drop = FALSE], 2, mu) }
  s <- apply(Xr[tr, , drop = FALSE], 2, sd); s[s == 0 | !is.finite(s)] <- 1; sweep(Xr, 2, s, "/") }
qtasks <- list(); for (r in 1:R) { f <- strat_folds(y, K, 1000 + r); for (k in 1:K) qtasks[[length(qtasks) + 1]] <- list(type = "random", r = r, k = k, te = which(f == k)) }
for (bb in levels(d$batch)) qtasks[[length(qtasks) + 1]] <- list(type = "lobo", b = bb, te = which(d$batch == bb))
qtasks[[length(qtasks) + 1]] <- list(type = "final", te = integer(0))
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261027)
clusterExport(cl, c("d", "X", "bres", "V")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
qres <- parLapply(cl, qtasks, function(t) { tr <- setdiff(seq_len(nrow(d)), t$te); Xr <- bres(X, d$batch, tr)
  m <- pipeline(d, Xr, tr); out <- t; out$genes <- m$genes
  if (length(t$te)) out$lp <- lp_new(m, Xr[t$te, , drop = FALSE]) else { out$lp61 <- lp_new(m, V$GSE63061$X); out$lp60 <- lp_new(m, V$GSE63060$X) }
  out })
stopCluster(cl)
oq <- matrix(NA, nrow(d), R); for (x in Filter(function(z) z$type == "random", qres)) oq[x$te, x$r] <- x$lp
lq <- rep(NA, nrow(d)); for (x in Filter(function(z) z$type == "lobo", qres)) lq[x$te] <- x$lp
strat <- function(s) { num <- 0; den <- 0; for (bb in levels(d$batch)) { i1 <- which(d$batch == bb & y == 1); i0 <- which(d$batch == bb & y == 0)
  if (length(i1) && length(i0)) { o <- outer(s[i1], s[i0], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) } }; num / den }
fq <- Filter(function(z) z$type == "final", qres)[[1]]
a_r <- sapply(1:R, function(r) auc1(y, oq[, r])); al <- auc_ci(y, lq); a61 <- auc_ci(V$GSE63061$d$AD, fq$lp61); a60 <- auc_ci(V$GSE63060$d$AD, fq$lp60)
Q <- data.frame(analysis = c("Random-fold CV, mean per-repeat AUC", "Random-fold CV, averaged-OOF batch-stratified AUC", "LOBO pooled AUC", "LOBO batch-stratified AUC",
                             "Final batch-residualised model: GSE63061 AUC", "Final batch-residualised model: GSE63060 AUC"),
  estimate = c(mean(a_r), strat(rowMeans(oq)), al[2], strat(lq), a61[2], a60[2]), lo = c(min(a_r), NA, al[1], NA, a61[1], a60[1]), hi = c(max(a_r), NA, al[3], NA, a61[3], a60[3]),
  note = c("lo/hi = range over 10 repeats", "", "DeLong", "", paste0(length(fq$genes), " genes; overlap with 17-gene model ", length(intersect(fq$genes, G0))), "DeLong"))
write.csv(Q, paste0(OUT, "Q_batch_residualised_development.csv"), row.names = FALSE); print(Q, digits = 3)

## ---- R. development uncertainty by subsampling without replacement (63.2%) ----
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261027)
clusterExport(cl, c("d", "X", "V")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
sub <- parLapply(cl, 1:500, function(b) { set.seed(9000 + b); si <- sample(seq_len(nrow(d)), round(0.632 * nrow(d)))
  m <- tryCatch(pipeline(d, X, si), error = function(e) NULL); if (is.null(m)) return(c(NA, NA, NA))
  c(as.numeric(auc(V$GSE63061$d$AD, lp_new(m, V$GSE63061$X), quiet = TRUE, direction = "<")),
    as.numeric(auc(V$GSE63060$d$AD, lp_new(m, V$GSE63060$X), quiet = TRUE, direction = "<")), length(m$genes)) })
stopCluster(cl)
sub <- do.call(rbind, sub); colnames(sub) <- c("AUC61", "AUC60", "n_genes"); write.csv(sub, paste0(OUT, "R_development_subsampling_draws.csv"), row.names = FALSE)
devb <- read.csv(paste0(OUT, "D_development_bootstrap_draws.csv"))
Rs <- do.call(rbind, lapply(list(bootstrap = devb, subsample_63.2pct = as.data.frame(sub)), function(z) data.frame(
  AUC61_median = median(z$AUC61, na.rm = TRUE), AUC61_lo = quantile(z$AUC61, 0.025, na.rm = TRUE), AUC61_hi = quantile(z$AUC61, 0.975, na.rm = TRUE),
  AUC60_median = median(z$AUC60, na.rm = TRUE), AUC60_lo = quantile(z$AUC60, 0.025, na.rm = TRUE), AUC60_hi = quantile(z$AUC60, 0.975, na.rm = TRUE),
  genes_median = median(z$n_genes, na.rm = TRUE), refits = sum(!is.na(z$AUC61)))))
Rs <- cbind(resampling = rownames(Rs), Rs); write.csv(Rs, paste0(OUT, "R_development_uncertainty_summary.csv"), row.names = FALSE); print(Rs, digits = 3)

## ---- S. random benchmarks from the step-A candidate pool and from strongly associated genes ----
top200 <- filter_genes(d, X, seq_len(nrow(d)), 200); poolA <- setdiff(top200, G0); poolS <- pool   # pool: |t| >= weakest signature gene (script 7)
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261027)
clusterExport(cl, c("d", "X", "V", "poolA", "poolS")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
bS <- parLapply(cl, 1:2000, function(i) { g <- sample(if (i <= 1000) poolA else poolS, 17); m <- fit_lrm(X, d$AD, seq_len(nrow(d)), g)
  c(type = ifelse(i <= 1000, 1, 2), a61 = as.numeric(auc(V$GSE63061$d$AD, lp_new(m, V$GSE63061$X), quiet = TRUE, direction = "<")),
    a60 = as.numeric(auc(V$GSE63060$d$AD, lp_new(m, V$GSE63060$X), quiet = TRUE, direction = "<"))) })
stopCluster(cl)
bS <- as.data.frame(do.call(rbind, bS)); bS$type <- ifelse(bS$type == 1, "step-A top-200 pool", "strong-association pool")
write.csv(bS, paste0(OUT, "S_random_benchmarks_candidate_pools_draws.csv"), row.names = FALSE)
pf <- read.csv("results/validation_performance.csv")
Ssum <- do.call(rbind, lapply(split(bS, bS$type), function(b) do.call(rbind, lapply(c("GSE63061", "GSE63060"), function(cn) {
  a <- if (cn == "GSE63061") b$a61 else b$a60; f <- pf$AUC[pf$cohort == cn]
  data.frame(benchmark = b$type[1], pool_size = ifelse(b$type[1] == "step-A top-200 pool", length(poolA), length(poolS)), cohort = cn,
             random_median = median(a), random_p95 = quantile(a, 0.95), model_AUC = f, percentile = mean(a < f) * 100, empirical_p = (sum(a >= f) + 1) / (length(a) + 1)) }))))
write.csv(Ssum, paste0(OUT, "S_random_benchmarks_candidate_pools_summary.csv"), row.names = FALSE); print(Ssum, digits = 3)

## ---- T. Brier score, scaled Brier, flexible (loess) calibration, ICI; age analyses ----
fv <- read.csv("results/validation_predictions.csv"); Tt <- list(); curves <- list(); ages <- list()
for (cn in names(V)) { v <- fv[fv$cohort == cn, ]; w <- V[[cn]]$d; v$age <- w$age[match(v$id, w$id)]
  br <- mean((v$p - v$AD)^2); br0 <- mean((mean(v$AD) - v$AD)^2)
  lo <- loess(AD ~ p, data = v, span = 0.75, degree = 2); pc <- pmin(pmax(predict(lo), 0), 1); ici <- mean(abs(pc - v$p))
  gr <- seq(quantile(v$p, 0.02), quantile(v$p, 0.98), length.out = 50)
  curves[[cn]] <- data.frame(cohort = cn, predicted = gr, observed_loess = pmin(pmax(predict(lo, data.frame(p = gr)), 0), 1))
  ## age: correlation of LP with age in controls / cases; AUC over age-matched pairs (|age difference| <= 2 y); age-residualised LP
  c0 <- cor.test(v$lp[v$AD == 0], v$age[v$AD == 0], method = "spearman", exact = FALSE); c1 <- cor.test(v$lp[v$AD == 1], v$age[v$AD == 1], method = "spearman", exact = FALSE)
  i1 <- which(v$AD == 1); i0 <- which(v$AD == 0); o <- outer(v$lp[i1], v$lp[i0], "-"); ad <- abs(outer(v$age[i1], v$age[i0], "-")) <= 2
  rl <- residuals(lm(lp ~ age, data = v[v$AD == 0, ])); fitc <- lm(lp ~ age, data = v[v$AD == 0, ]); lr <- v$lp - predict(fitc, v)
  Tt[[cn]] <- data.frame(cohort = cn, Brier = br, Brier_null = br0, scaled_Brier = 1 - br / br0, ICI_loess = ici,
    rho_LP_age_controls = unname(c0$estimate), p_rho_controls = c0$p.value, rho_LP_age_AD = unname(c1$estimate), p_rho_AD = c1$p.value,
    AUC_age_matched_pairs = (sum(o[ad] > 0) + 0.5 * sum(o[ad] == 0)) / sum(ad), n_age_matched_pairs = sum(ad),
    AUC_age_residualised = auc1(v$AD, lr), AUC_age_only = auc1(v$AD, v$age)) }
lpo <- read.csv("results/discovery_oof_predictions.csv"); lpd <- lpo$oof_lp_mean[match(d$id, lpo$id)]
c0 <- cor.test(lpd[y == 0], d$age[y == 0], method = "spearman", exact = FALSE)
Tt[["GSE140829"]] <- data.frame(cohort = "GSE140829 (averaged OOF)", Brier = mean((plogis(lpd) - y)^2), Brier_null = mean((mean(y) - y)^2),
  scaled_Brier = 1 - mean((plogis(lpd) - y)^2) / mean((mean(y) - y)^2), ICI_loess = NA, rho_LP_age_controls = unname(c0$estimate), p_rho_controls = c0$p.value,
  rho_LP_age_AD = NA, p_rho_AD = NA, AUC_age_matched_pairs = NA, n_age_matched_pairs = NA, AUC_age_residualised = NA, AUC_age_only = auc1(y, d$age))
Tt <- do.call(rbind, Tt); write.csv(Tt, paste0(OUT, "T_brier_flexible_calibration_age.csv"), row.names = FALSE); print(Tt, digits = 3)
write.csv(do.call(rbind, curves), paste0(OUT, "T_loess_calibration_curves.csv"), row.names = FALSE)

## ---- U. sex-discordant samples retained (AddNeuroMed) ----
U <- list(); bdraw <- list()
for (cn in names(V)) { w <- read_cohort(sprintf("data/%s_inclSexDiscordant_z.csv.gz", cn)); Xw <- as.matrix(w[, genes]); lp <- lp_new(fm$final, Xw)
  a <- auc_ci(w$AD, lp); cc <- cal(w$AD, lp); dc <- w$sex_discordant %in% c(TRUE, "True", "TRUE")
  set.seed(31); ra <- replicate(1000, { g <- character(0); for (p in ex_pool) g <- c(g, sample(setdiff(p, g), 1)); m <- fit_lrm(X, y, seq_len(nrow(d)), g); auc1(w$AD, lp_new(m, Xw)) })
  wt <- if (sum(dc & w$AD == 0) > 1) wilcox.test(lp[dc & w$AD == 0], lp[!dc & w$AD == 0]) else NULL
  U[[cn]] <- data.frame(cohort = cn, n = nrow(w), n_AD = sum(w$AD), n_CN = sum(1 - w$AD), n_sex_discordant = sum(dc), discordant_CN = sum(dc & w$AD == 0),
    AUC = a[2], lo = a[1], hi = a[3], cal_intercept = cc["intercept"], cal_slope = cc["slope"],
    random_expression_matched_percentile = mean(ra < a[2]) * 100, random_median = median(ra),
    mean_LP_discordant_CN = mean(lp[dc & w$AD == 0]), mean_LP_concordant_CN = mean(lp[!dc & w$AD == 0]), mean_LP_AD = mean(lp[w$AD == 1]),
    wilcoxon_p_discordant_vs_concordant_CN = if (is.null(wt)) NA else wt$p.value,
    AUC_concordant_CN_only_vs_all_AD = auc1(c(w$AD[w$AD == 1], w$AD[!dc & w$AD == 0]), c(lp[w$AD == 1], lp[!dc & w$AD == 0]))) }
U <- do.call(rbind, U); write.csv(U, paste0(OUT, "U_sex_discordant_retained.csv"), row.names = FALSE); print(U, digits = 3)

## ---- V. lambda-rule use across outer folds; per-batch AUC random-effects summary (LOBO) ----
sl <- read.csv("results/discovery_cv_fold_selection_log.csv"); rl <- unique(sl[, c("rep", "fold", "rule", "penalty")])
Vr <- as.data.frame(table(rule = rl$rule)); Vr$share <- Vr$Freq / sum(Vr$Freq)
write.csv(Vr, paste0(OUT, "V_lambda_rule_counts_outer_folds.csv"), row.names = FALSE); print(Vr)
pb <- read.csv(paste0(OUT, "C_lobo_per_batch.csv")); pe <- pb[!is.na(pb$AUC), ]
hm <- with(pe, { q1 <- AUC / (2 - AUC); q2 <- 2 * AUC^2 / (1 + AUC); (AUC * (1 - AUC) + (n_AD - 1) * (q1 - AUC^2) + (n_CN - 1) * (q2 - AUC^2)) / (n_AD * n_CN) })
Lg <- qlogis(pe$AUC); vl <- hm / (pe$AUC * (1 - pe$AUC))^2; wf <- 1 / vl; fe <- sum(wf * Lg) / sum(wf); Qc <- sum(wf * (Lg - fe)^2); dfq <- nrow(pe) - 1   # logit-AUC scale
tau2 <- max(0, (Qc - dfq) / (sum(wf) - sum(wf^2) / sum(wf))); wr <- 1 / (vl + tau2); re <- sum(wr * Lg) / sum(wr); se <- sqrt(1 / sum(wr)); tq <- qt(0.975, max(1, dfq - 1))
Vb <- data.frame(batches = nrow(pe), RE_AUC = plogis(re), lo = plogis(re - 1.96 * se), hi = plogis(re + 1.96 * se), tau2_logit = tau2, I2 = max(0, (Qc - dfq) / Qc), Q = Qc, Q_p = pchisq(Qc, dfq, lower.tail = FALSE),
  pred_lo = plogis(re - tq * sqrt(se^2 + tau2)), pred_hi = plogis(re + tq * sqrt(se^2 + tau2)))
write.csv(Vb, paste0(OUT, "V_lobo_per_batch_random_effects.csv"), row.names = FALSE); print(Vb, digits = 3)

## ---- W. composition-only model developed in discovery and applied unchanged; EPIC deconvolution sensitivity ----
CSd <- as.data.frame(cscore(X)); dd <- cbind(d[, c("AD", "age", "sex")], CSd)
mc <- glm(AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK, binomial, dd)
dd$lp_sig <- lpd; mcs <- glm(AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK + lp_sig, binomial, dd)
mc0 <- glm(AD ~ Neutrophil + Monocyte + Bcell + Tcell + NK, binomial, dd); mcs0 <- glm(AD ~ Neutrophil + Monocyte + Bcell + Tcell + NK + lp_sig, binomial, dd)   # without age/sex (age effect differs between cohorts)
write.csv(data.frame(term = names(coef(mc)), coefficient = coef(mc)), paste0(OUT, "W_composition_model_coefficients.csv"), row.names = FALSE)
W <- list()
for (cn in names(V)) { w <- cbind(V[[cn]]$d[, c("AD", "age", "sex")], as.data.frame(cscore(V[[cn]]$X))); w$lp_sig <- lp_new(fm$final, V[[cn]]$X)
  w$lp_comp <- predict(mc, w); w$lp_both <- predict(mcs, w); w$lp_cell <- predict(mc0, w); w$lp_cellsig <- predict(mcs0, w)
  g0 <- glm(AD ~ lp_comp, binomial, w); g1 <- glm(AD ~ lp_comp + lp_sig, binomial, w); lr <- anova(g0, g1, test = "LRT")
  ac <- auc_ci(w$AD, w$lp_comp); ab <- auc_ci(w$AD, w$lp_both); as_ <- auc_ci(w$AD, w$lp_sig); a0c <- auc_ci(w$AD, w$lp_cell); a0b <- auc_ci(w$AD, w$lp_cellsig)
  h0 <- glm(AD ~ lp_cell, binomial, w); h1 <- glm(AD ~ lp_cell + lp_sig, binomial, w); lr0 <- anova(h0, h1, test = "LRT")
  W[[cn]] <- data.frame(cohort = cn, AUC_signature = as_[2], AUC_composition_model = ac[2], comp_lo = ac[1], comp_hi = ac[3],
    AUC_composition_plus_signature = ab[2], both_lo = ab[1], both_hi = ab[3], LR_chisq_signature_beyond_composition = lr$Deviance[2], LR_p = lr$`Pr(>Chi)`[2],
    cal_slope_composition = cal(w$AD, w$lp_comp)["slope"], cor_sig_comp = cor(w$lp_sig, w$lp_comp),
    AUC_cells_only = a0c[2], cells_lo = a0c[1], cells_hi = a0c[3], AUC_cells_plus_signature = a0b[2], cellsig_lo = a0b[1], cellsig_hi = a0b[3],
    LR_chisq_signature_beyond_cells = lr0$Deviance[2], LR_p_cells = lr0$`Pr(>Chi)`[2], cor_sig_cells = cor(w$lp_sig, w$lp_cell)) }
W <- do.call(rbind, W); write.csv(W, paste0(OUT, "W_frozen_composition_model.csv"), row.names = FALSE); print(W, digits = 3)
.libPaths(c(normalizePath("../Rlib", mustWork = FALSE), normalizePath("Rlib", mustWork = FALSE), .libPaths()))
if (requireNamespace("EPIC", quietly = TRUE)) {
  E <- list()
  for (cn in c("GSE140829", "GSE63061", "GSE63060")) {
    L <- read.csv(gzfile(sprintf("data/%s_log2_allgenes.csv.gz", cn)), check.names = FALSE); ids <- L$id; M <- t(2^as.matrix(L[, -1])); colnames(M) <- ids
    ep <- EPIC::EPIC(bulk = M, reference = EPIC::BRef, withOtherCells = TRUE)$cellFractions
    w <- if (cn == "GSE140829") d else V[[cn]]$d; fr <- as.data.frame(ep[match(w$id, rownames(ep)), c("Bcells", "CD4_Tcells", "CD8_Tcells", "Monocytes", "Neutrophils", "NKcells")])
    write.csv(cbind(id = w$id, fr), sprintf("%sW_EPIC_fractions_%s.csv", OUT, cn), row.names = FALSE)
    w <- cbind(w[, intersect(c("AD", "age", "sex", "batch"), names(w))], fr); w$lpz <- as.numeric(scale(if (cn == "GSE140829") lpd else lp_new(fm$final, V[[cn]]$X)))
    f0 <- AD ~ age + sex + Bcells + CD4_Tcells + CD8_Tcells + Monocytes + Neutrophils + NKcells; if (cn == "GSE140829") f0 <- update(f0, . ~ . + batch)
    m0 <- glm(f0, binomial, w); m1 <- update(m0, . ~ . + lpz); b <- coef(summary(m1))["lpz", ]
    E[[cn]] <- data.frame(cohort = cn, OR_adjusted_EPIC = exp(b[1]), lo = exp(b[1] - 1.96 * b[2]), hi = exp(b[1] + 1.96 * b[2]), LR_p = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2],
      AUC_covariates = auc1(w$AD, fitted(m0)), AUC_cov_plus_sig = auc1(w$AD, fitted(m1)),
      mean_neutrophil_fraction = mean(fr$Neutrophils), EPIC_genes_used = sum(EPIC::BRef$sigGenes %in% colnames(L))) }
  E <- do.call(rbind, E); write.csv(E, paste0(OUT, "W_adjusted_EPIC_deconvolution.csv"), row.names = FALSE); print(E, digits = 3) }

## ---- X. specificity contrasts: delta AUC (AD vs CN minus group vs CN), bootstrap CI, minimum detectable difference ----
ref <- "Control (GSE140829)"; s_ref <- scM[grpAll == ref]; s_ad <- scM[grpAll == "AD (GSE140829)"]; set.seed(41)
hmse <- function(A, n1, n0) { q1 <- A / (2 - A); q2 <- 2 * A^2 / (1 + A); sqrt((A * (1 - A) + (n1 - 1) * (q1 - A^2) + (n0 - 1) * (q2 - A^2)) / (n1 * n0)) }
Xs <- do.call(rbind, lapply(gr2, function(g) { s_g <- scM[grpAll == g]
  a_ad <- pairauc(s_ad, s_ref); a_g <- pairauc(s_g, s_ref)
  bt <- replicate(2000, { r0 <- sample(s_ref, replace = TRUE); pairauc(sample(s_ad, replace = TRUE), r0) - pairauc(sample(s_g, replace = TRUE), r0) })
  se_g <- hmse(0.5, length(s_g), length(s_ref))
  data.frame(group = g, n = length(s_g), AUC_AD_vs_CN = a_ad, AUC_group_vs_CN = a_g, delta = a_ad - a_g, delta_lo = quantile(bt, 0.025), delta_hi = quantile(bt, 0.975),
    MDD_AUC_vs_0.5 = 2.80 * se_g, MDD_note = "smallest AUC departure from 0.5 detectable with 80% power (two-sided alpha 0.05), Hanley-McNeil SE") }))
write.csv(Xs, paste0(OUT, "X_specificity_delta_auc.csv"), row.names = FALSE); print(Xs, digits = 3)

## ---- Y. APOE e4 and the signature (discovery, averaged OOF score) ----
if ("apoe" %in% names(d)) { e4 <- ifelse(is.na(d$apoe) | d$apoe == "", NA, as.integer(grepl("4", d$apoe)))
  dd2 <- data.frame(AD = y, lpz = as.numeric(scale(lpd)), e4 = e4, age = d$age, sex = d$sex, batch = d$batch)
  m_a <- glm(AD ~ lpz + age + sex + batch, binomial, dd2[!is.na(e4), ]); m_b <- update(m_a, . ~ . + e4); m_c <- update(m_b, . ~ . + lpz:e4)
  t_cn <- t.test(dd2$lpz[y == 0 & e4 %in% 1], dd2$lpz[y == 0 & e4 %in% 0]); t_ad <- t.test(dd2$lpz[y == 1 & e4 %in% 1], dd2$lpz[y == 1 & e4 %in% 0])
  Y <- data.frame(metric = c("APOE genotype available (n)", "e4 carriers AD (%)", "e4 carriers CN (%)", "OR/SD signature, age+sex+batch", "OR/SD signature, + e4",
    "LR p signature x e4 interaction", "Mean score difference e4+ vs e4- in CN (SD)", "p", "Mean score difference e4+ vs e4- in AD (SD)", "p", "AUC e4 alone", "AUC averaged-OOF score (same subset)"),
    value = c(sum(!is.na(e4)), 100 * mean(e4[y == 1], na.rm = TRUE), 100 * mean(e4[y == 0], na.rm = TRUE), exp(coef(m_a)["lpz"]), exp(coef(m_b)["lpz"]),
      anova(m_b, m_c, test = "LRT")$`Pr(>Chi)`[2], diff(rev(t_cn$estimate)), t_cn$p.value, diff(rev(t_ad$estimate)), t_ad$p.value,
      auc1(y[!is.na(e4)], e4[!is.na(e4)]), auc1(y[!is.na(e4)], lpd[!is.na(e4)])))
  write.csv(Y, paste0(OUT, "Y_apoe.csv"), row.names = FALSE); print(Y, digits = 3) }

## ---- Z. probe-level transfer: identical-sequence probes and probe averaging ----
Z <- list()
for (cn in names(V)) for (tg in c("fixedprobe", "probemean")) { f <- sprintf("data/%s_%s_z.csv.gz", cn, tg); if (!file.exists(f)) next
  z <- read.csv(gzfile(f), check.names = FALSE); w <- V[[cn]]$d; z <- z[match(w$id, z$id), ]; Xz <- as.matrix(z[, G0]); colnames(Xz) <- G0
  lp <- lp_new(fm$final, Xz); a <- auc_ci(w$AD, lp); cc <- cal(w$AD, lp)
  Z[[paste(cn, tg)]] <- data.frame(cohort = cn, probes = ifelse(tg == "fixedprobe", "probe with identical sequence to discovery probe", "mean of all probes annotated to gene"),
    AUC = a[2], lo = a[1], hi = a[3], cal_intercept = cc["intercept"], cal_slope = cc["slope"], cor_with_primary_LP = cor(lp, lp_new(fm$final, V[[cn]]$X))) }
Z <- do.call(rbind, Z); write.csv(Z, paste0(OUT, "Z_probe_level_transfer.csv"), row.names = FALSE); print(Z, digits = 3)
writeLines(c(sprintf("discovery prevalence,%.4f", mean(y)), sprintf("threshold probability,%.4f", fm$threshold)), paste0(OUT, "Z_threshold_vs_prevalence.csv"))
writeLines(capture.output(sessionInfo()), paste0(OUT, "sessionInfo_R.txt"))
