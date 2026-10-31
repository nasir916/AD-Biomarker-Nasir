
## 7_secondary_analyses.R — further secondary analyses (uncertainty, specificity, benchmarks, supporting statistics)
## (run after 5_sensitivity_analyses.R in the same session, or it re-sources the needed objects)
if (!exists("res")) source("5_sensitivity_analyses.R")
suppressPackageStartupMessages(library(parallel)); set.seed(20261026)
OUT <- "results/sensitivity/"; y <- DAT$cohortz$disc$AD
fv <- read.csv("results/validation_predictions.csv")
## ---- 1. non-stratified participant bootstrap + model-based CIs for validation metrics ----
nb <- list()
for (cn in c("GSE63061", "GSE63060")) {
  key <- ifelse(cn == "GSE63061", "v61", "v60"); v <- fv[fv$cohort == cn, ]; vv <- DAT$cohortz[[key]]
  v$age <- vv$age[match(v$id, vv$id)]; v$sex <- vv$sex[match(v$id, vv$id)]
  st <- function(w) { a0 <- glm(AD ~ age + sex, binomial, w); a1 <- glm(AD ~ age + sex + lp, binomial, w)
    c(as.numeric(auc(w$AD, w$lp, quiet = TRUE, direction = "<")), cal(w$AD, w$lp),
      as.numeric(auc(w$AD, fitted(a1), quiet = TRUE, direction = "<")) - as.numeric(auc(w$AD, fitted(a0), quiet = TRUE, direction = "<"))) }
  S <- t(replicate(2000, st(v[sample(nrow(v), replace = TRUE), ])))
  g0 <- glm(AD ~ offset(lp), binomial, v); g1 <- glm(AD ~ lp, binomial, v)
  wi <- confint.default(g0)[1, ]; ws <- confint.default(g1)[2, ]
  est <- st(v)
  nb[[cn]] <- data.frame(cohort = cn, metric = c("AUC", "Calibration intercept", "Calibration slope", "Delta AUC (apparent, in-cohort)"),
    estimate = est, boot_lo = apply(S, 2, quantile, 0.025), boot_hi = apply(S, 2, quantile, 0.975),
    wald_lo = c(NA, wi[1], ws[1], NA), wald_hi = c(NA, wi[2], ws[2], NA)) }
nb <- do.call(rbind, nb); write.csv(nb, paste0(OUT, "I_validation_uncertainty_nonstratified.csv"), row.names = FALSE); print(nb, digits = 3)
## ---- 2. threshold on the final model's scale ----
fm <- readRDS("results/frozen_model.rds"); lp_app <- lp_new(fm$final, XM$cohortz$disc); p_app <- plogis(lp_app)
oofd <- read.csv("results/discovery_oof_predictions.csv"); p_oof <- plogis(oofd$oof_lp_mean)
thr_t <- data.frame(metric = c("threshold", "fraction of discovery above threshold: averaged OOF", "fraction above threshold: final model (apparent)",
  "apparent sensitivity (final model)", "apparent specificity (final model)", "SD of averaged OOF LP", "SD of final-model apparent LP", "mean per-repeat OOF LP SD"),
  value = c(fm$threshold, mean(p_oof >= fm$threshold), mean(p_app >= fm$threshold), mean(p_app[y == 1] >= fm$threshold), mean(p_app[y == 0] < fm$threshold),
            sd(oofd$oof_lp_mean), sd(lp_app), NA))
write.csv(thr_t, paste0(OUT, "J_threshold_scale.csv"), row.names = FALSE); print(thr_t, digits = 3)
## ---- 3. specificity: matched ensembles (one fold model per repeat for every participant), AD vs MCI, batch-stratified ----
rr <- Filter(function(x) x$cn == "fixedref" && x$r > 0, res)
oofF <- matrix(NA, length(y), R); for (x in rr) oofF[x$te, x$r] <- x$lp_te
set.seed(11); extM <- matrix(NA, nrow(EXT), R)
for (r in 1:R) { pick <- sample(1:K, nrow(EXT), replace = TRUE)
  for (k in 1:K) { x <- Filter(function(z) z$r == r && z$k == k, rr)[[1]]; ii <- which(pick == k); extM[ii, r] <- x$lp_ext[ii] } }
grpAll <- c(ifelse(y == 1, "AD (GSE140829)", "Control (GSE140829)"), EXT$grp); scM <- c(rowMeans(oofF), rowMeans(extM))
## batch labels for every scored participant
b_disc <- as.character(DAT$cohortz$disc$batch)
mci_b <- read.csv(gzfile("data/GSE140829_MCI_fixedref.csv.gz"), check.names = FALSE)[, c("id", "batch")]
ftd_raw <- read.csv(gzfile("data/GSE140830_fixedref.csv.gz"), check.names = FALSE)
dupids <- read.csv("results/GSE140830_duplicates_of_GSE140829.csv"); ftd_raw <- ftd_raw[!(ftd_raw$id %in% dupids$id[dupids$dup_of != ""]), ]
b_ext <- c(if ("batch" %in% names(ftd_raw)) as.character(ftd_raw$batch) else rep(NA, nrow(ftd_raw)), as.character(mci_b$batch))
batchAll <- c(b_disc, b_ext)
cmp2 <- function(a, b) { s1 <- scM[grpAll == a]; s0 <- scM[grpAll == b]; ci <- auc_ci(c(rep(1, length(s1)), rep(0, length(s0))), c(s1, s0))
  ## batch-stratified AUC (pairs within the same batch only)
  i1 <- which(grpAll == a); i0 <- which(grpAll == b); num <- 0; den <- 0
  for (bb in intersect(unique(batchAll[i1]), unique(batchAll[i0]))) { if (is.na(bb)) next
    o <- outer(scM[i1[batchAll[i1] == bb]], scM[i0[batchAll[i0] == bb]], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) }
  data.frame(case = a, reference = b, n_case = length(s1), n_ref = length(s0), AUC = ci[2], lo = ci[1], hi = ci[3],
             batch_stratified_AUC = if (den > 0) num / den else NA, within_batch_pairs = den) }
gr2 <- c("bvFTD", "PSP", "nfvPPA", "svPPA", "CBS", "MCI (GSE140829)")
sp2 <- rbind(cmp2("AD (GSE140829)", "Control (GSE140829)"), do.call(rbind, lapply(gr2, function(g) cmp2(g, "Control (GSE140829)"))),
             do.call(rbind, lapply(gr2, function(g) cmp2("AD (GSE140829)", g))), cmp2("Independent control (GSE140830)", "Control (GSE140829)"))
write.csv(sp2, paste0(OUT, "K_specificity_matched_ensemble.csv"), row.names = FALSE); print(sp2, digits = 3)
bt <- table(group = grpAll, batch = batchAll, useNA = "ifany"); write.csv(as.data.frame.matrix(bt), paste0(OUT, "K_specificity_batch_overlap.csv"))
## frozen 17-gene model, fixed-reference scaling, external groups vs the 33 independent controls (no training data involved)
lpF <- lp_new(fm$final, XE); s0 <- lpF[EXT$grp == "Independent control (GSE140830)"]
fz <- do.call(rbind, lapply(gr2, function(g) { s1 <- lpF[EXT$grp == g]; ci <- auc_ci(c(rep(1, length(s1)), rep(0, length(s0))), c(s1, s0))
  data.frame(case = g, reference = paste0("Independent control (GSE140830), n = ", length(s0)), n_case = length(s1), AUC = ci[2], lo = ci[1], hi = ci[3]) }))
write.csv(fz, paste0(OUT, "K_specificity_frozen_model_vs_independent_controls.csv"), row.names = FALSE); print(fz, digits = 3)
## ---- 4. random-signature benchmarks with exactly 17 genes: expression-matched and composition-matched ----
Xd <- XM$cohortz$disc
mk <- list(Neutrophil = c("FCGR3B","CSF3R","CXCR2","FPR1"), Monocyte = c("CD14","LYZ","CSF1R","FCN1"), Bcell = c("MS4A1","CD19","CD79A","BANK1"),
           Tcell = c("CD3E","CD3D","CD2","LCK"), NK = c("GNLY","NKG7","KLRD1","PRF1"))
CS <- sapply(mk, function(g) rowMeans(Xd[, intersect(g, colnames(Xd))]))
gcor <- cor(Xd, CS); mu <- read.csv("results/discovery_gene_means.csv"); dec <- setNames(mu$bin, mu$gene)
G0 <- fm$genes; cand <- setdiff(colnames(Xd), G0)
nn_pool <- lapply(G0, function(g) { pool <- cand[abs(dec[cand] - dec[g]) <= 1]
  d2 <- rowSums((gcor[pool, , drop = FALSE] - matrix(gcor[g, ], length(pool), 5, byrow = TRUE))^2); pool[order(d2)][1:50] })
ex_pool <- lapply(G0, function(g) cand[dec[cand] == dec[g]])
draw17 <- function(pools) { out <- character(0); for (p in pools) { p <- setdiff(p, out); out <- c(out, sample(p, 1)) }; out }
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261026)
clusterExport(cl, c("DAT", "XM", "nn_pool", "ex_pool", "draw17")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
bench2 <- parLapply(cl, 1:2000, function(i) { pools <- if (i <= 1000) ex_pool else nn_pool; g <- draw17(pools)
  m <- fit_lrm(XM$cohortz$disc, DAT$cohortz$disc$AD, seq_len(nrow(DAT$cohortz$disc)), g)
  c(type = ifelse(i <= 1000, 1, 2), n = length(unique(g)),
    a61 = as.numeric(auc(DAT$cohortz$v61$AD, lp_new(m, XM$cohortz$v61), quiet = TRUE, direction = "<")),
    a60 = as.numeric(auc(DAT$cohortz$v60$AD, lp_new(m, XM$cohortz$v60), quiet = TRUE, direction = "<"))) })
stopCluster(cl)
bench2 <- as.data.frame(do.call(rbind, bench2)); bench2$type <- ifelse(bench2$type == 1, "expression-matched", "composition-matched")
write.csv(bench2, paste0(OUT, "L_random_benchmarks_17genes_draws.csv"), row.names = FALSE)
pf <- read.csv("results/validation_performance.csv")
bsum <- do.call(rbind, lapply(split(bench2, bench2$type), function(b) do.call(rbind, lapply(c("GSE63061", "GSE63060"), function(cn) {
  a <- if (cn == "GSE63061") b$a61 else b$a60; f <- pf$AUC[pf$cohort == cn]
  data.frame(benchmark = b$type[1], cohort = cn, genes_per_set = paste(range(b$n), collapse = "-"), random_median = median(a), random_p95 = quantile(a, 0.95),
             frozen_AUC = f, percentile = mean(a < f) * 100, empirical_p = (sum(a >= f) + 1) / (length(a) + 1)) }))))
write.csv(bsum, paste0(OUT, "L_random_benchmarks_17genes_summary.csv"), row.names = FALSE); print(bsum, digits = 3)
## ---- 5. direction-concordance null: random 17-gene sets drawn from genes with discovery |t| at least as large as the weakest signature gene ----
tt <- read.csv("results/adjusted_tstats_all_cohorts.csv"); rownames(tt) <- tt$gene
thr_t <- min(abs(tt[G0, "t_GSE140829"])); pool <- setdiff(tt$gene[abs(tt$t_GSE140829) >= thr_t], G0)
conc <- function(gs, col) sum(sign(tt[gs, "t_GSE140829"]) == sign(tt[gs, col]))
nullc <- t(replicate(5000, { g <- sample(pool, 17); c(conc(g, "t_GSE63061"), conc(g, "t_GSE63060")) }))
cn_out <- data.frame(cohort = c("GSE63061", "GSE63060"), pool_size = length(pool), discovery_abs_t_min = thr_t,
  observed = c(conc(G0, "t_GSE63061"), conc(G0, "t_GSE63060")), null_mean = colMeans(nullc),
  empirical_p = c((sum(nullc[, 1] >= conc(G0, "t_GSE63061")) + 1) / 5001, (sum(nullc[, 2] >= conc(G0, "t_GSE63060")) + 1) / 5001))
write.csv(cn_out, paste0(OUT, "M_concordance_null.csv"), row.names = FALSE); print(cn_out, digits = 3)
## ---- 6. supporting statistics ----
d <- DAT$cohortz$disc; ct <- table(d$batch, d$AD); chi <- suppressWarnings(chisq.test(ct)); V <- sqrt(chi$statistic / (sum(ct) * (min(dim(ct)) - 1)))
mk6 <- c(mk, list(Platelet_erythroid = c("PPBP", "PF4", "HBB", "ALAS2")))
CS6 <- sapply(mk6, function(g) rowMeans(Xd[, intersect(g, colnames(Xd)), drop = FALSE]))
eff <- do.call(rbind, lapply(colnames(CS6), function(k) { m <- lm(CS6[, k] ~ AD + age + sex + batch, data = d); s <- summary(m)$coefficients["AD", ]
  data.frame(score = k, genes_used = paste(intersect(mk6[[k]], colnames(Xd)), collapse = ","), AD_beta = s[1], p = s[4]) }))
gc6 <- round(cor(Xd[, G0], CS6), 3)
write.csv(gc6, paste0(OUT, "N_signature_gene_celltype_correlations.csv")); write.csv(eff, paste0(OUT, "N_celltype_score_AD_effects.csv"), row.names = FALSE)
lpo <- oofd$oof_lp_mean[match(d$id, oofd$id)]; num <- 0; den <- 0
for (b in levels(d$batch)) { i1 <- which(d$batch == b & y == 1); i0 <- which(d$batch == b & y == 0)
  if (length(i1) && length(i0)) { o <- outer(lpo[i1], lpo[i0], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) } }
per_rep_strat <- sapply(1:R, function(r) { lr <- NULL; NA })
sup <- data.frame(metric = c("Batch-diagnosis Cramer's V (discovery)", "Batch-diagnosis chi-square p", "Random-fold OOF (averaged LP) batch-stratified AUC",
  "R2 of averaged OOF LP on batch", "Mean per-repeat random-fold CV AUC", "SD of per-repeat CV AUC", "Discovery age p (Student t)", "Discovery age p (Welch t)"),
  value = c(V, chi$p.value, num / den, summary(lm(lpo ~ d$batch))$r.squared, mean(read.csv("results/discovery_cv_per_repeat.csv")$auc_signature),
    sd(read.csv("results/discovery_cv_per_repeat.csv")$auc_signature), t.test(d$age[y == 1], d$age[y == 0], var.equal = TRUE)$p.value, t.test(d$age[y == 1], d$age[y == 0])$p.value))
write.csv(sup, paste0(OUT, "N_supporting_statistics.csv"), row.names = FALSE); print(sup, digits = 3); print(eff, digits = 3)
## ---- 7. alternative marker sets for confounding adjustment ----
mk2 <- list(Neutrophil = c("CXCR1","MMP9","S100A12","CEACAM8","FCGR3B","CSF3R","CXCR2","FPR1"), Monocyte = c("VCAN","CD163","FCGR1A","CD14","LYZ","CSF1R","FCN1"),
            Bcell = c("CD22","PAX5","FCRL5","MS4A1","CD19","BANK1"), Tcell = c("CD3G","CD5","CD6","CD28","CD3E","CD3D","CD2","LCK"), NK = c("KLRF1","NCAM1","SH2D1B","GZMB","GNLY","NKG7","KLRD1","PRF1"))
adj2 <- list()
for (cn in c("GSE140829", "GSE63061", "GSE63060")) {
  key <- switch(cn, GSE140829 = "disc", GSE63061 = "v61", GSE63060 = "v60"); w <- DAT$cohortz[[key]]; Xk <- XM$cohortz[[key]]
  for (k in names(mk2)) w[[k]] <- rowMeans(Xk[, intersect(setdiff(mk2[[k]], G0), colnames(Xk)), drop = FALSE])
  w$lp <- if (cn == "GSE140829") oofd$oof_lp_mean[match(w$id, oofd$id)] else fv$lp[fv$cohort == cn][match(w$id, fv$id[fv$cohort == cn])]
  w$lpz <- as.numeric(scale(w$lp)); f0 <- AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK
  if (cn == "GSE140829") f0 <- update(f0, . ~ . + batch)
  m0 <- glm(f0, binomial, w); m1 <- update(m0, . ~ . + lpz); b <- coef(summary(m1))["lpz", ]
  adj2[[cn]] <- data.frame(cohort = cn, markers = "extended (6-8 genes per cell type, signature genes excluded)", OR_adjusted = exp(b[1]),
    lo = exp(b[1] - 1.96 * b[2]), hi = exp(b[1] + 1.96 * b[2]), LR_p = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2],
    genes_used = paste(sapply(mk2, function(g) length(intersect(setdiff(g, G0), colnames(Xk)))), collapse = "/")) }
adj2 <- do.call(rbind, adj2); write.csv(adj2, paste0(OUT, "O_adjusted_alternative_markers.csv"), row.names = FALSE); print(adj2, digits = 3)
