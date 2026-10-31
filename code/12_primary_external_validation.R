
## 12_primary_external_validation.R — primary external analysis (fixed-reference scaling; sex-discordant samples retained), its
## random-signature benchmarks, site/sex/age/composition analyses, age- and sex-adjusted specificity, LOBO with the held-out batch
## standardised on itself, and an unweighted-score comparator. Run after 5_, 7_ and 10_ in the same R session (run_all.py).
if (!exists("scM") || !exists("poolA")) { source("5_sensitivity_analyses.R"); source("7_secondary_analyses.R"); source("10_additional_sensitivity_analyses.R") }
suppressPackageStartupMessages(library(parallel)); OUT <- "results/sensitivity/"; set.seed(20261028)
fm <- readRDS("results/frozen_model.rds"); G0 <- fm$genes; thr <- fm$threshold
d <- DAT$cohortz$disc; X <- XM$cohortz$disc; y <- d$AD
auc1 <- function(yy, s) as.numeric(auc(yy, s, quiet = TRUE, direction = "<"))
hmse <- function(A, n1, n0) { q1 <- A / (2 - A); q2 <- 2 * A^2 / (1 + A); sqrt((A * (1 - A) + (n1 - 1) * (q1 - A^2) + (n0 - 1) * (q2 - A^2)) / (n1 * n0)) }
pairauc <- function(s1, s0) { o <- outer(s1, s0, "-"); mean(o > 0) + 0.5 * mean(o == 0) }
re_auc <- function(A, n1, n0) {   # DerSimonian-Laird random-effects summary on the logit-AUC scale (Hanley-McNeil SE, delta method)
  A <- pmin(pmax(A, 0.005), 0.995); L <- qlogis(A); v <- (hmse(A, n1, n0) / (A * (1 - A)))^2; w <- 1 / v; fe <- sum(w * L) / sum(w); Qc <- sum(w * (L - fe)^2); k <- length(A)
  tau2 <- max(0, (Qc - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); wr <- 1 / (v + tau2); re <- sum(wr * L) / sum(wr); se <- sqrt(1 / sum(wr)); tq <- qt(0.975, max(1, k - 2))
  c(k = k, RE_AUC = plogis(re), lo = plogis(re - 1.96 * se), hi = plogis(re + 1.96 * se), tau2_logit = tau2, I2 = max(0, (Qc - (k - 1)) / Qc),
    pred_lo = plogis(re - tq * sqrt(se^2 + tau2)), pred_hi = plogis(re + tq * sqrt(se^2 + tau2))) }
## ---------- load validation sets ----------
ld <- function(cn, tag) { x <- read_cohort(sprintf("data/%s_%s.csv.gz", cn, tag)); x$AD <- as.integer(x$dx %in% c("AD")); x }
SETS <- list()
for (cn in c("GSE63061", "GSE63060")) {
  inf <- read.csv(gzfile(sprintf("data/%s_inclSexDiscordant_z.csv.gz", cn)), check.names = FALSE)[, c("id", "sex_inferred", "sex_discordant")]
  for (tag in c("retained_fixedref", "fixedref", "retained_cohortz", "cohortz", "retained_ctlref", "ctlref")) {
    w <- ld(cn, tag); w$sex_inf <- factor(inf$sex_inferred[match(w$id, inf$id)], levels = c("female", "male"))
    w$discordant <- inf$sex_discordant[match(w$id, inf$id)] %in% c(TRUE, "True", "TRUE")
    SETS[[paste(cn, tag)]] <- w } }
Xof <- function(w) as.matrix(w[, genes])
prim <- list(GSE63061 = SETS[["GSE63061 retained_fixedref"]], GSE63060 = SETS[["GSE63060 retained_fixedref"]])
## ---------- A. primary performance ----------
m2 <- read.csv("results/frozen_model2_lp_age_sex.csv"); b2 <- setNames(m2$coefficient, m2$term)
as0 <- glm(AD ~ age + sex, binomial, d)        # age + sex model developed in discovery
PA <- list(); PRED <- list(); CUR <- list()
for (cn in names(prim)) { w <- prim[[cn]]; w$lp <- lp_new(fm$final, Xof(w)); w$p <- plogis(w$lp); PRED[[cn]] <- data.frame(cohort = cn, id = w$id, AD = w$AD, lp = w$lp, p = w$p)
  a <- auc_ci(w$AD, w$lp); g0 <- glm(AD ~ offset(lp), binomial, w); g1 <- glm(AD ~ lp, binomial, w); wi <- confint.default(g0)[1, ]; ws <- confint.default(g1)[2, ]
  br <- mean((w$p - w$AD)^2); br0 <- mean((mean(w$AD) - w$AD)^2); lo <- loess(AD ~ p, data = w, span = 0.75, degree = 2)
  gr <- seq(quantile(w$p, 0.02), quantile(w$p, 0.98), length.out = 50); CUR[[cn]] <- data.frame(cohort = cn, predicted = gr, observed_loess = pmin(pmax(predict(lo, data.frame(p = gr)), 0), 1))
  pc <- as.integer(w$p >= thr); TP <- sum(pc == 1 & w$AD == 1); FN <- sum(pc == 0 & w$AD == 1); FP <- sum(pc == 1 & w$AD == 0); TN <- sum(pc == 0 & w$AD == 0)
  ws1 <- wilson(TP, TP + FN); ws0 <- wilson(TN, TN + FP)
  w$sex <- w$sex_inf; a0 <- glm(AD ~ age + sex, binomial, w); a1 <- glm(AD ~ age + sex + lp, binomial, w); lr <- anova(a0, a1, test = "LRT")
  dA <- function(ww) { f0 <- glm(AD ~ age + sex, binomial, ww); f1 <- glm(AD ~ age + sex + lp, binomial, ww); auc1(ww$AD, fitted(f1)) - auc1(ww$AD, fitted(f0)) }
  bt <- replicate(2000, dA(w[sample(nrow(w), replace = TRUE), ]))
  lp_as <- predict(as0, newdata = w); lp_c <- b2["(Intercept)"] + b2["lp"] * w$lp + b2["age"] * w$age + b2[grep("sex", names(b2))] * (w$sex == "male")
  PA[[cn]] <- data.frame(cohort = cn, n = nrow(w), n_AD = sum(w$AD), n_CN = sum(1 - w$AD), prevalence = mean(w$AD), AUC = a[2], AUC_lo = a[1], AUC_hi = a[3],
    CITL = coef(g0)[1], CITL_lo = wi[1], CITL_hi = wi[2], slope = coef(g1)[2], slope_lo = ws[1], slope_hi = ws[2], slope_CI_width = ws[2] - ws[1],
    Brier = br, Brier_null = br0, scaled_Brier = 1 - br / br0, ICI = mean(abs(pmin(pmax(predict(lo), 0), 1) - w$p)),
    TP = TP, FN = FN, FP = FP, TN = TN, sens = TP / (TP + FN), sens_lo = ws1[1], sens_hi = ws1[2], spec = TN / (TN + FP), spec_lo = ws0[1], spec_hi = ws0[2],
    PPV = TP / (TP + FP), NPV = TN / (TN + FN), LR_chisq_over_age_sex = lr$Deviance[2], LR_p = lr$`Pr(>Chi)`[2],
    AUC_incohort_age_sex = auc1(w$AD, fitted(a0)), AUC_incohort_age_sex_sig = auc1(w$AD, fitted(a1)), dAUC = dA(w), dAUC_lo = quantile(bt, 0.025), dAUC_hi = quantile(bt, 0.975),
    AUC_age_sex_model_from_discovery = auc1(w$AD, lp_as), AUC_combined_model_from_discovery = auc1(w$AD, lp_c), row.names = NULL) }
PA <- do.call(rbind, PA); write.csv(PA, paste0(OUT, "AA_primary_external_performance.csv"), row.names = FALSE); print(t(PA), digits = 3)
write.csv(do.call(rbind, PRED), paste0(OUT, "AA_primary_external_predictions.csv"), row.names = FALSE); write.csv(do.call(rbind, CUR), paste0(OUT, "AA_primary_loess_curves.csv"), row.names = FALSE)
## ---------- B. random-signature benchmarks (exactly 17 genes) on the primary external data ----------
XP <- lapply(prim, Xof); YP <- lapply(prim, function(w) w$AD)
draw17 <- function(pools) { out <- character(0); for (p in pools) { p <- setdiff(p, out); out <- c(out, sample(p, 1)) }; out }
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261028)
clusterExport(cl, c("d", "X", "XP", "YP", "ex_pool", "nn_pool", "poolA", "poolS", "draw17")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
BB <- parLapply(cl, 1:4000, function(i) { ty <- (i - 1) %/% 1000 + 1
  g <- switch(ty, draw17(ex_pool), draw17(nn_pool), sample(poolA, 17), sample(poolS, 17))
  m <- fit_lrm(X, d$AD, seq_len(nrow(d)), g)
  c(type = ty, a61 = as.numeric(auc(YP$GSE63061, lp_new(m, XP$GSE63061), quiet = TRUE, direction = "<")), a60 = as.numeric(auc(YP$GSE63060, lp_new(m, XP$GSE63060), quiet = TRUE, direction = "<"))) })
stopCluster(cl)
BB <- as.data.frame(do.call(rbind, BB)); BB$type <- c("expression-matched", "composition-matched", "step-A candidate pool", "strong-association pool")[BB$type]
write.csv(BB, paste0(OUT, "AB_primary_random_benchmarks_draws.csv"), row.names = FALSE)
BS <- do.call(rbind, lapply(split(BB, BB$type), function(b) do.call(rbind, lapply(names(prim), function(cn) { a <- if (cn == "GSE63061") b$a61 else b$a60; f <- PA$AUC[PA$cohort == cn]
  data.frame(benchmark = b$type[1], cohort = cn, random_median = median(a), random_p95 = quantile(a, 0.95), model_AUC = f, percentile = mean(a < f) * 100, empirical_p = (sum(a >= f) + 1) / (length(a) + 1)) }))))
write.csv(BS, paste0(OUT, "AB_primary_random_benchmarks_summary.csv"), row.names = FALSE); print(BS, digits = 3)
## ---------- C. scaling x analysis-set grid ----------
CG <- do.call(rbind, lapply(names(SETS), function(k) { w <- SETS[[k]]; lp <- lp_new(fm$final, Xof(w)); a <- auc_ci(w$AD, lp); cc <- cal(w$AD, lp); p <- plogis(lp)
  data.frame(cohort = strsplit(k, " ")[[1]][1], set = ifelse(grepl("retained", k), "sex-discordant retained", "sex-discordant excluded"),
    scaling = sub("retained_", "", strsplit(k, " ")[[1]][2]), n_AD = sum(w$AD), n_CN = sum(1 - w$AD), AUC = a[2], lo = a[1], hi = a[3], CITL = cc["intercept"], slope = cc["slope"],
    Brier = mean((p - w$AD)^2), scaled_Brier = 1 - mean((p - w$AD)^2) / mean((mean(w$AD) - w$AD)^2)) }))
write.csv(CG, paste0(OUT, "AC_scaling_by_analysis_set.csv"), row.names = FALSE); print(CG, digits = 3)
## ---------- D. sites (primary data) ----------
site_of <- function(cn) { md <- readLines(gzfile(sprintf("raw/%s.gz", cn)), n = 400)
  acc <- strsplit(md[grep("^!Sample_geo_accession", md)], "\t")[[1]][-1]; des <- strsplit(md[grep("^!Sample_description", md)][1], "\t")[[1]][-1]
  setNames(toupper(substr(gsub('"', "", des), 1, 3)), gsub('"', "", acc)) }
SD <- list(); SR <- list()
for (cn in names(prim)) { w <- prim[[cn]]; w$lp <- lp_new(fm$final, Xof(w)); w$site <- site_of(cn)[w$id]
  ps <- do.call(rbind, lapply(split(w, w$site), function(s) data.frame(cohort = cn, site = s$site[1], n_AD = sum(s$AD), n_CN = sum(1 - s$AD),
    AUC = if (sum(s$AD) >= 5 && sum(1 - s$AD) >= 5) auc1(s$AD, s$lp) else NA)))
  SD[[cn]] <- ps; e <- ps[!is.na(ps$AUC), ]; r <- re_auc(e$AUC, e$n_AD, e$n_CN)
  num <- 0; den <- 0; for (s in unique(w$site)) { i1 <- which(w$site == s & w$AD == 1); i0 <- which(w$site == s & w$AD == 0); if (length(i1) && length(i0)) { o <- outer(w$lp[i1], w$lp[i0], "-"); num <- num + sum(o > 0) + 0.5 * sum(o == 0); den <- den + length(o) } }
  SR[[cn]] <- data.frame(cohort = cn, t(r), site_stratified_AUC = num / den) }
write.csv(do.call(rbind, SD), paste0(OUT, "AD_primary_per_site_AUC.csv"), row.names = FALSE); SR <- do.call(rbind, SR); write.csv(SR, paste0(OUT, "AD_primary_site_random_effects.csv"), row.names = FALSE); print(SR, digits = 3)
## ---------- E/F. sex-stratified and age analyses (primary data; GSE63060 control metadata may be unreliable) ----------
oofd <- read.csv("results/discovery_oof_predictions.csv"); lpd <- oofd$oof_lp_mean[match(d$id, oofd$id)]
EF <- list()
for (cn in names(prim)) { w <- prim[[cn]]; w$lp <- lp_new(fm$final, Xof(w))
  i1 <- which(w$AD == 1); i0 <- which(w$AD == 0); o <- outer(w$lp[i1], w$lp[i0], "-"); am <- abs(outer(w$age[i1], w$age[i0], "-")) <= 2
  fc <- lm(lp ~ age, data = w[w$AD == 0, ]); c0 <- cor.test(w$lp[i0], w$age[i0], method = "spearman", exact = FALSE)
  EF[[cn]] <- data.frame(cohort = cn, AUC_female = auc1(w$AD[w$sex_inf == "female"], w$lp[w$sex_inf == "female"]), n_female = sum(w$sex_inf == "female"),
    AUC_male = auc1(w$AD[w$sex_inf == "male"], w$lp[w$sex_inf == "male"]), n_male = sum(w$sex_inf == "male"),
    rho_age_controls = unname(c0$estimate), p_rho = c0$p.value, AUC_age_matched_pairs = (sum(o[am] > 0) + 0.5 * sum(o[am] == 0)) / sum(am), n_pairs = sum(am),
    AUC_age_residualised = auc1(w$AD, w$lp - predict(fc, w)), age_AD = mean(w$age[i1]), age_CN = mean(w$age[i0]), p_age = t.test(w$age[i1], w$age[i0], var.equal = TRUE)$p.value) }
EF[["GSE140829 (OOF)"]] <- data.frame(cohort = "GSE140829 (averaged OOF)", AUC_female = auc1(y[d$sex == "female"], lpd[d$sex == "female"]), n_female = sum(d$sex == "female"),
  AUC_male = auc1(y[d$sex == "male"], lpd[d$sex == "male"]), n_male = sum(d$sex == "male"), rho_age_controls = cor(lpd[y == 0], d$age[y == 0], method = "spearman"),
  p_rho = cor.test(lpd[y == 0], d$age[y == 0], method = "spearman", exact = FALSE)$p.value, AUC_age_matched_pairs = NA, n_pairs = NA, AUC_age_residualised = NA,
  age_AD = mean(d$age[y == 1]), age_CN = mean(d$age[y == 0]), p_age = t.test(d$age[y == 1], d$age[y == 0], var.equal = TRUE)$p.value)
EF <- do.call(rbind, EF); write.csv(EF, paste0(OUT, "AE_primary_sex_age_analyses.csv"), row.names = FALSE); print(EF, digits = 3)
## ---------- G. cell composition on the primary data ----------
GC <- list()
for (cn in names(prim)) { w <- prim[[cn]]; Xw <- Xof(w); w$lp <- lp_new(fm$final, Xw); w$sex <- w$sex_inf
  cs <- as.data.frame(cscore(Xw)); w <- cbind(w[, c("AD", "age", "sex", "lp")], cs); w$lpz <- as.numeric(scale(w$lp))
  m0 <- glm(AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK, binomial, w); m1 <- update(m0, . ~ . + lpz); b <- coef(summary(m1))["lpz", ]
  w$lp_cell <- predict(mc0, w); w$lp_comp <- predict(mc, w); w$lp_sig <- w$lp
  h0 <- glm(AD ~ lp_cell, binomial, w); h1 <- update(h0, . ~ . + lp_sig); k0 <- glm(AD ~ lp_comp, binomial, w); k1 <- update(k0, . ~ . + lp_sig)
  w2 <- w; for (k in names(mk2)) w2[[k]] <- rowMeans(Xw[, intersect(setdiff(mk2[[k]], G0), colnames(Xw)), drop = FALSE])
  e0 <- glm(AD ~ age + sex + Neutrophil + Monocyte + Bcell + Tcell + NK, binomial, w2); e1 <- update(e0, . ~ . + lpz); be <- coef(summary(e1))["lpz", ]
  L <- read.csv(gzfile(sprintf("data/%s_retained_log2_allgenes.csv.gz", cn)), check.names = FALSE); M <- t(2^as.matrix(L[, -1])); colnames(M) <- L$id
  ep <- EPIC::EPIC(bulk = M, reference = EPIC::BRef, withOtherCells = TRUE)$cellFractions; fr <- as.data.frame(ep[match(prim[[cn]]$id, rownames(ep)), c("Bcells", "CD4_Tcells", "CD8_Tcells", "Monocytes", "Neutrophils", "NKcells")])
  w3 <- cbind(w[, c("AD", "age", "sex", "lpz")], fr); f0 <- glm(AD ~ age + sex + Bcells + CD4_Tcells + CD8_Tcells + Monocytes + Neutrophils + NKcells, binomial, w3); f1 <- update(f0, . ~ . + lpz); bf <- coef(summary(f1))["lpz", ]
  GC[[cn]] <- data.frame(cohort = cn, OR_markers = exp(b[1]), OR_markers_lo = exp(b[1] - 1.96 * b[2]), OR_markers_hi = exp(b[1] + 1.96 * b[2]), LR_p_markers = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2],
    AUC_markers_incohort = auc1(w$AD, fitted(m0)), AUC_markers_plus_sig_incohort = auc1(w$AD, fitted(m1)),
    OR_extended = exp(be[1]), OR_extended_lo = exp(be[1] - 1.96 * be[2]), OR_extended_hi = exp(be[1] + 1.96 * be[2]), LR_p_extended = anova(e0, e1, test = "LRT")$`Pr(>Chi)`[2],
    OR_EPIC = exp(bf[1]), OR_EPIC_lo = exp(bf[1] - 1.96 * bf[2]), OR_EPIC_hi = exp(bf[1] + 1.96 * bf[2]), LR_p_EPIC = anova(f0, f1, test = "LRT")$`Pr(>Chi)`[2],
    AUC_cells_model_from_discovery = auc1(w$AD, w$lp_cell), AUC_cells_lo = auc_ci(w$AD, w$lp_cell)[1], AUC_cells_hi = auc_ci(w$AD, w$lp_cell)[3], LR_p_sig_beyond_cells = anova(h0, h1, test = "LRT")$`Pr(>Chi)`[2],
    AUC_comp_model_from_discovery = auc1(w$AD, w$lp_comp), LR_p_sig_beyond_comp = anova(k0, k1, test = "LRT")$`Pr(>Chi)`[2],
    AUC_cells_plus_sig_from_discovery = auc1(w$AD, predict(mcs0, w)), AUC_comp_plus_sig_from_discovery = auc1(w$AD, predict(mcs, w))) }
GC <- do.call(rbind, GC); write.csv(GC, paste0(OUT, "AG_primary_composition.csv"), row.names = FALSE); print(t(GC), digits = 3)
## ---------- H. specificity adjusted for age and sex ----------
ftd_m <- read_cohort("data/GSE140830_fixedref.csv.gz"); dupx <- read.csv("results/GSE140830_duplicates_of_GSE140829.csv"); ftd_m <- ftd_m[!(ftd_m$id %in% dupx$id[dupx$dup_of != ""]), ]
mci_m <- read_cohort("data/GSE140829_MCI_fixedref.csv.gz")
ageAll <- c(d$age, ftd_m$age, mci_m$age); sexAll <- factor(tolower(c(as.character(d$sex), as.character(ftd_m$sex), as.character(mci_m$sex))), levels = c("female", "male"))
stopifnot(length(ageAll) == length(scM)); ref <- "Control (GSE140829)"; i_ref <- which(grpAll == ref)
fc <- lm(s ~ age + sex, data = data.frame(s = scM[i_ref], age = ageAll[i_ref], sex = sexAll[i_ref]))      # age/sex dependence estimated in discovery controls
scR <- scM - predict(fc, data.frame(age = ageAll, sex = sexAll))
GT <- do.call(rbind, lapply(unique(grpAll), function(g) { i <- grpAll == g; data.frame(group = g, n = sum(i), age_mean = mean(ageAll[i], na.rm = TRUE), age_sd = sd(ageAll[i], na.rm = TRUE),
  age_missing = sum(is.na(ageAll[i])), female_pct = 100 * mean(sexAll[i] == "female", na.rm = TRUE)) }))
write.csv(GT, paste0(OUT, "AH_specificity_group_age_sex.csv"), row.names = FALSE); print(GT, digits = 3)
grs <- c("AD (GSE140829)", "PSP", "nfvPPA", "CBS", "bvFTD", "svPPA", "MCI (GSE140829)", "Independent control (GSE140830)")
s_ad_r <- scR[grpAll == "AD (GSE140829)"]; s_ref_r <- scR[i_ref]
SA <- do.call(rbind, lapply(grs, function(g) { i1 <- which(grpAll == g); dd <- data.frame(case = c(rep(1, length(i1)), rep(0, length(i_ref))), s = as.numeric(scale(c(scM[i1], scM[i_ref]))),
    age = ageAll[c(i1, i_ref)], sex = sexAll[c(i1, i_ref)])
  m0 <- glm(case ~ age + sex, binomial, dd); m1 <- update(m0, . ~ . + s); b <- coef(summary(m1))["s", ]
  o <- outer(scM[i1], scM[i_ref], "-"); am <- abs(outer(ageAll[i1], ageAll[i_ref], "-")) <= 2
  ar <- auc_ci(c(rep(1, length(i1)), rep(0, length(i_ref))), c(scR[i1], scR[i_ref]))
  bt <- if (g == "AD (GSE140829)") c(NA, NA) else quantile(replicate(2000, { r0 <- sample(s_ref_r, replace = TRUE); pairauc(sample(s_ad_r, replace = TRUE), r0) - pairauc(sample(scR[i1], replace = TRUE), r0) }), c(0.025, 0.975))
  data.frame(group = g, n = length(i1), AUC_unadjusted = pairauc(scM[i1], scM[i_ref]), AUC_age_sex_residualised = ar[2], resid_lo = ar[1], resid_hi = ar[3],
    AUC_age_matched_pairs = (sum(o[am] > 0) + 0.5 * sum(o[am] == 0)) / sum(am), n_age_matched_pairs = sum(am),
    OR_per_SD_age_sex_adjusted = exp(b[1]), OR_lo = exp(b[1] - 1.96 * b[2]), OR_hi = exp(b[1] + 1.96 * b[2]), LR_p = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2],
    dAUC_residualised = if (g == "AD (GSE140829)") NA else pairauc(s_ad_r, s_ref_r) - pairauc(scR[i1], s_ref_r), dAUC_lo = bt[1], dAUC_hi = bt[2]) }))
SA$LR_p_BH <- NA; jj <- SA$group != "AD (GSE140829)"; SA$LR_p_BH[jj] <- p.adjust(SA$LR_p[jj], "BH")
write.csv(SA, paste0(OUT, "AH_specificity_age_sex_adjusted.csv"), row.names = FALSE); print(SA, digits = 3)
## direct comparisons AD vs each syndrome on the age/sex-residualised score
SD2 <- do.call(rbind, lapply(grs[2:7], function(g) { s1 <- scR[grpAll == "AD (GSE140829)"]; s0 <- scR[grpAll == g]; ci <- auc_ci(c(rep(1, length(s1)), rep(0, length(s0))), c(s1, s0))
  data.frame(comparison = paste("AD vs", g), AUC_residualised = ci[2], lo = ci[1], hi = ci[3]) }))
write.csv(SD2, paste0(OUT, "AH_specificity_AD_vs_group_residualised.csv"), row.names = FALSE); print(SD2, digits = 3)
## unchanged 17-gene model vs the independent controls, with age/sex adjustment
lpF <- lp_new(fm$final, XE); iC <- which(EXT$grp == "Independent control (GSE140830)"); ageE <- c(ftd_m$age, mci_m$age); sexE <- sexAll[-(1:nrow(d))]
FZ <- do.call(rbind, lapply(c("PSP", "nfvPPA", "CBS", "bvFTD", "svPPA", "MCI (GSE140829)"), function(g) { i1 <- which(EXT$grp == g)
  dd <- data.frame(case = c(rep(1, length(i1)), rep(0, length(iC))), s = as.numeric(scale(lpF[c(i1, iC)])), age = ageE[c(i1, iC)], sex = sexE[c(i1, iC)])
  m0 <- glm(case ~ age + sex, binomial, dd); m1 <- update(m0, . ~ . + s); b <- coef(summary(m1))["s", ]
  data.frame(group = g, n = length(i1), n_ref = length(iC), AUC = pairauc(lpF[i1], lpF[iC]), OR_per_SD_age_sex_adjusted = exp(b[1]), OR_lo = exp(b[1] - 1.96 * b[2]), OR_hi = exp(b[1] + 1.96 * b[2]), LR_p = anova(m0, m1, test = "LRT")$`Pr(>Chi)`[2]) }))
write.csv(FZ, paste0(OUT, "AH_specificity_unchanged_model_vs_independent_controls_adjusted.csv"), row.names = FALSE); print(FZ, digits = 3)
## FDR within the cell-type family
ce <- read.csv(paste0(OUT, "N_celltype_score_AD_effects.csv")); ce$p_BH <- p.adjust(ce$p, "BH"); write.csv(ce, paste0(OUT, "AH_celltype_effects_BH.csv"), row.names = FALSE); print(ce)
## ---------- I. LOBO with the held-out batch standardised on itself (training batches as in the main analysis) ----------
cl <- makeCluster(NCORE); clusterSetRNGStream(cl, 20261028); clusterExport(cl, c("d", "X")); clusterEvalQ(cl, { source("R/pipeline_functions.R"); NULL })
LB <- parLapply(cl, levels(d$batch), function(b) { tr <- which(d$batch != b); te <- which(d$batch == b); m <- pipeline(d, X, tr)
  Xt <- X[te, , drop = FALSE]; if (length(te) > 2) { Xt <- scale(Xt); Xt[is.na(Xt)] <- 0 }; list(b = b, te = te, lp = lp_new(m, Xt)) })
stopCluster(cl)
lps <- rep(NA, nrow(d)); for (x in LB) lps[x$te] <- x$lp; pl <- auc_ci(y, lps)
pb2 <- do.call(rbind, lapply(LB, function(x) { yy <- y[x$te]; data.frame(batch = x$b, n_AD = sum(yy), n_CN = sum(1 - yy), AUC = if (sum(yy) >= 5 && sum(1 - yy) >= 5) auc1(yy, x$lp) else NA) }))
e <- pb2[!is.na(pb2$AUC), ]; rr <- re_auc(e$AUC, e$n_AD, e$n_CN)
LBO <- data.frame(analysis = "LOBO, held-out batch standardised on itself", pooled_AUC = pl[2], lo = pl[1], hi = pl[3], batch_stratified_AUC = strat(lps), t(rr))
pbm <- read.csv(paste0(OUT, "C_lobo_per_batch.csv")); pbm <- pbm[!is.na(pbm$AUC), ]; rm_ <- re_auc(pbm$AUC, pbm$n_AD, pbm$n_CN)
LBO <- rbind(LBO, data.frame(analysis = "LOBO, main analysis (per-batch AUCs)", pooled_AUC = NA, lo = NA, hi = NA, batch_stratified_AUC = NA, t(rm_)))
write.csv(LBO, paste0(OUT, "AI_lobo_heldout_self_standardised.csv"), row.names = FALSE); print(LBO, digits = 3)
## ---------- J. unweighted-score comparators ----------
fa <- filter_genes(d, X, seq_len(nrow(d)), 200); tt0 <- read.csv("results/adjusted_tstats_all_cohorts.csv"); rownames(tt0) <- tt0$gene
top17 <- fa[1:17]; sg <- sign(tt0[top17, "t_GSE140829"]); sg0 <- sign(fm$final$fit$coefficients[make.names(G0)])
UJ <- do.call(rbind, lapply(names(prim), function(cn) { Xw <- XP[[cn]]
  u1 <- as.numeric(Xw[, top17] %*% sg) / 17; u2 <- as.numeric(Xw[, G0] %*% sign(tt0[G0, "t_GSE140829"])) / 17; a1 <- auc_ci(YP[[cn]], u1); a2 <- auc_ci(YP[[cn]], u2)
  data.frame(cohort = cn, AUC_unweighted_top17_stepA = a1[2], lo1 = a1[1], hi1 = a1[3], AUC_unweighted_signature_genes = a2[2], lo2 = a2[1], hi2 = a2[3], overlap_top17_with_signature = length(intersect(top17, G0))) }))
write.csv(UJ, paste0(OUT, "AJ_unweighted_scores.csv"), row.names = FALSE); write.csv(data.frame(gene = top17, sign = sg), paste0(OUT, "AJ_top17_stepA_genes.csv"), row.names = FALSE); print(UJ, digits = 3)
writeLines(capture.output(sessionInfo()), paste0(OUT, "sessionInfo_R.txt"))
