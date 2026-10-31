
## 2_validation.R — frozen-model external validation in GSE63061 (primary) and GSE63060 (sensitivity)
## Methods 2.6–2.7. No refitting of the signature in validation data.
suppressPackageStartupMessages({library(rms); library(pROC); library(parallel)})
set.seed(20261025)
fm <- readRDS("results/frozen_model.rds"); final <- fm$final; thr <- fm$threshold; G <- fm$genes
meta_cols <- c("id","AD","age","sex","batch","apoe")
rd <- function(f) { x <- read.csv(gzfile(f), check.names = FALSE); x$sex <- factor(tolower(x$sex), levels = c("female","male")); x }
disc <- rd("data/GSE140829_z.csv.gz")
cohorts <- list(GSE63061 = rd("data/GSE63061_z.csv.gz"), GSE63060 = rd("data/GSE63060_z.csv.gz"))
universe <- setdiff(names(disc), meta_cols)
lp_of <- function(m, df, genes) {
  nd <- as.data.frame(df[, genes, drop = FALSE]); names(nd) <- make.names(genes)
  for (g in setdiff(make.names(genes), names(nd))) nd[[g]] <- 0
  predict(m, newdata = nd, type = "lp") }
## secondary frozen model: signature LP (out-of-fold in discovery) + age + sex
oofd <- read.csv("results/discovery_oof_predictions.csv")
disc$lp <- oofd$oof_lp_mean[match(disc$id, oofd$id)]
m2 <- glm(AD ~ lp + age + sex, binomial, disc); m_as <- glm(AD ~ age + sex, binomial, disc)
write.csv(data.frame(term = names(coef(m2)), coefficient = coef(m2)), "results/frozen_model2_lp_age_sex.csv", row.names = FALSE)

perf <- list(); preds <- list()
for (nm in names(cohorts)) {
  v <- cohorts[[nm]]; y <- v$AD
  missing <- setdiff(G, names(v))
  v$lp <- lp_of(final$fit, v, G); p <- plogis(v$lp)
  rc <- roc(y, v$lp, quiet = TRUE, direction = "<"); ci <- ci.auc(rc, method = "delong")
  cal_slope <- coef(glm(y ~ v$lp, binomial))[2]
  cal_citl <- coef(glm(y ~ offset(v$lp), binomial))[1]
  pred_cls <- as.integer(p >= thr)
  tp <- sum(pred_cls == 1 & y == 1); tn <- sum(pred_cls == 0 & y == 0); fp <- sum(pred_cls == 1 & y == 0); fn <- sum(pred_cls == 0 & y == 1)
  # frozen secondary model and frozen age+sex model
  p2 <- predict(m2, v, type = "link"); pas <- predict(m_as, v, type = "link")
  rc2 <- roc(y, p2, quiet = TRUE, direction = "<"); rcas <- roc(y, pas, quiet = TRUE, direction = "<")
  # added value in-cohort: LR test of LP over age+sex
  a0 <- glm(y ~ age + sex, binomial, v); a1 <- glm(y ~ age + sex + lp, binomial, v)
  lr <- anova(a0, a1, test = "LRT")
  r0 <- roc(y, fitted(a0), quiet = TRUE, direction = "<"); r1 <- roc(y, fitted(a1), quiet = TRUE, direction = "<")
  dl <- roc.test(r1, r0, method = "delong", paired = TRUE)
  perf[[nm]] <- data.frame(cohort = nm, n = length(y), n_AD = sum(y), n_CN = sum(1 - y), genes_missing = length(missing),
    AUC = ci[2], AUC_lo = ci[1], AUC_hi = ci[3], cal_intercept_citl = cal_citl, cal_slope = cal_slope,
    sens = tp / (tp + fn), spec = tn / (tn + fp), PPV = tp / (tp + fp), NPV = tn / (tn + fn),
    AUC_frozen_age_sex = as.numeric(auc(rcas)), AUC_frozen_model2 = as.numeric(ci.auc(rc2)[2]),
    AUC_frozen_model2_lo = ci.auc(rc2)[1], AUC_frozen_model2_hi = ci.auc(rc2)[3],
    LR_chisq_LP_over_age_sex = lr$Deviance[2], LR_p = lr$`Pr(>Chi)`[2],
    AUC_incohort_age_sex = as.numeric(auc(r0)), AUC_incohort_age_sex_LP = as.numeric(auc(r1)), dAUC_p = dl$p.value,
    TP = tp, TN = tn, FP = fp, FN = fn)
  preds[[nm]] <- data.frame(cohort = nm, id = v$id, AD = y, lp = v$lp, p = p, p_model2 = plogis(p2))
}
perf <- do.call(rbind, perf); rownames(perf) <- NULL
write.csv(perf, "results/validation_performance.csv", row.names = FALSE)
write.csv(do.call(rbind, preds), "results/validation_predictions.csv", row.names = FALSE)
print(t(perf))

## agreement: age/sex-adjusted AD t-statistics, genome-wide and for signature genes
tstat <- function(df) {
  D <- model.matrix(~ AD + age + sex, df); Y <- as.matrix(df[, universe]); q <- qr(D)
  B <- qr.coef(q, Y); Rr <- Y - D %*% B; s2 <- colSums(Rr^2) / (nrow(D) - q$rank)
  XtXi <- chol2inv(qr.R(q)); B[2, ] / sqrt(s2 * XtXi[2, 2]) }
t_d <- tstat(disc); agree <- list(); tt <- data.frame(gene = universe, t_GSE140829 = t_d)
sgn_coef <- sign(coef(final$fit)[-1]); names(sgn_coef) <- G
for (nm in names(cohorts)) {
  t_v <- tstat(cohorts[[nm]]); tt[[paste0("t_", nm)]] <- t_v
  rho <- cor.test(t_d, t_v, method = "spearman", exact = FALSE)
  conc <- sum(sign(t_d[G]) == sign(t_v[G])); bt <- binom.test(conc, length(G), 0.5, alternative = "greater")
  concc <- sum(sgn_coef == sign(t_v[G])); btc <- binom.test(concc, length(G), 0.5, alternative = "greater")
  agree[[nm]] <- data.frame(cohort = nm, spearman_rho_genomewide = rho$estimate, rho_p = rho$p.value,
    sig_genes = length(G), concordant_marginal = conc, sign_test_p = bt$p.value,
    concordant_with_model_coef = concc, sign_test_coef_p = btc$p.value)
}
agree <- do.call(rbind, agree); write.csv(agree, "results/validation_agreement.csv", row.names = FALSE)
write.csv(tt, "results/adjusted_tstats_all_cohorts.csv", row.names = FALSE); print(agree)

## random-signature benchmarks (exactly 17 genes) are computed in 7_secondary_analyses.R and 12_primary_external_validation.R
