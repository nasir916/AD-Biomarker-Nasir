
## 1_discovery.R — signature development in GSE140829 (Methods 2.4–2.5)
suppressPackageStartupMessages({library(glmnet); library(rms); library(pROC)})
set.seed(20261025)
d <- read.csv(gzfile("data/GSE140829_z.csv.gz"), check.names = FALSE)
meta_cols <- c("id","AD","age","sex","batch","apoe")
genes_all <- setdiff(names(d), meta_cols)
d$sex <- factor(tolower(d$sex)); d$batch <- factor(d$batch)
X_all <- as.matrix(d[, genes_all]); y <- d$AD

filter_genes <- function(idx, top = 200) {
  D <- model.matrix(~ AD + age + sex + batch, data = droplevels(d[idx, ]))
  Y <- X_all[idx, , drop = FALSE]
  q <- qr(D); B <- qr.coef(q, Y); R <- Y - D %*% B
  df <- nrow(D) - q$rank; s2 <- colSums(R^2) / df
  XtXi <- chol2inv(qr.R(q))[1:q$rank, 1:q$rank]
  j <- which(colnames(D) == "AD")
  t <- B[j, ] / sqrt(s2 * XtXi[j, j])
  p <- 2 * pt(-abs(t), df)
  names(sort(p))[1:top]
}
select_en <- function(idx, g) {
  cv <- cv.glmnet(X_all[idx, g], y[idx], family = "binomial", alpha = 0.5, nfolds = 10)
  nz <- function(l) { b <- coef(cv, s = l)[-1, 1]; names(b)[b != 0] }
  sel <- nz(cv$lambda.1se); rule <- "1se"
  if (length(sel) > 20) {               # least-penalised lambda with <= 20 non-zero genes
    ok <- which(cv$glmnet.fit$df <= 20); l <- cv$glmnet.fit$lambda[max(ok)]
    sel <- nz(l); rule <- "df<=20" }
  if (length(sel) == 0) { sel <- nz(cv$lambda.min); rule <- "min" }
  list(genes = sel, rule = rule, cv = cv)
}
fit_lrm <- function(idx, g) {
  dd <- data.frame(AD = y[idx], X_all[idx, g, drop = FALSE], check.names = FALSE)
  names(dd) <- make.names(names(dd))
  f <- as.formula(paste("AD ~", paste(make.names(g), collapse = " + ")))
  f0 <- lrm(f, data = dd, x = TRUE, y = TRUE)
  pt <- pentrace(f0, c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128))
  pen <- pt$penalty
  fit <- if (pen > 0) update(f0, penalty = pen) else f0
  list(fit = fit, penalty = pen, genes = g)
}
lp_new <- function(m, newX) {
  nd <- as.data.frame(newX[, m$genes, drop = FALSE]); names(nd) <- make.names(m$genes)
  predict(m$fit, newdata = nd, type = "lp")
}
pipeline <- function(idx) { g <- filter_genes(idx); s <- select_en(idx, g); m <- fit_lrm(idx, s$genes); m$rule <- s$rule; m$cv <- s$cv; m$candidates <- g; m }

## --- repeated nested CV (5 outer folds x 10 repeats) ---
R <- 10; K <- 5; oof <- matrix(NA, nrow(d), R); sel_log <- list(); fold_auc <- c(); FOLDS <- matrix(NA, nrow(d), R); fold_slope <- c(); fold_citl <- c()
oof_as <- oof_las <- matrix(NA, nrow(d), R)
for (r in 1:R) {
  folds <- rep(NA, nrow(d))
  for (cl in 0:1) { ii <- which(y == cl); folds[ii] <- sample(rep(1:K, length.out = length(ii))) }
  FOLDS[, r] <- folds
  for (k in 1:K) {
    tr <- which(folds != k); te <- which(folds == k)
    m <- pipeline(tr)
    lp_tr <- lp_new(m, X_all[tr, ]); oof[te, r] <- lp_new(m, X_all[te, ])
    sel_log[[length(sel_log) + 1]] <- data.frame(rep = r, fold = k, gene = m$genes, rule = m$rule, penalty = m$penalty)
    # reference age+sex and LP+age+sex fitted on training fold only
    dtr <- data.frame(AD = y[tr], age = d$age[tr], sex = d$sex[tr], lp = lp_tr)
    dte <- data.frame(age = d$age[te], sex = d$sex[te], lp = oof[te, r])
    a1 <- glm(AD ~ age + sex, binomial, dtr); a2 <- glm(AD ~ lp + age + sex, binomial, dtr)
    oof_as[te, r] <- predict(a1, dte); oof_las[te, r] <- predict(a2, dte)
    fold_auc <- c(fold_auc, as.numeric(auc(y[te], oof[te, r], quiet = TRUE)))
    fold_slope <- c(fold_slope, unname(coef(glm(y[te] ~ oof[te, r], binomial))[2]))          # fold model evaluated on its own held-out fold
    fold_citl <- c(fold_citl, unname(coef(glm(y[te] ~ offset(oof[te, r]), binomial))[1]))
  }
  cat("repeat", r, "AUC", round(as.numeric(auc(y, oof[, r], quiet = TRUE)), 3), "\n")
}
sel <- do.call(rbind, sel_log)
per_rep <- function(M) sapply(1:R, function(r) as.numeric(auc(y, M[, r], quiet = TRUE)))
slope <- sapply(1:R, function(r) coef(glm(y ~ oof[, r], binomial))[2])
citl  <- sapply(1:R, function(r) coef(glm(y ~ offset(oof[, r]), binomial))[1])
cv_summary <- data.frame(
  metric = c("CV-AUC signature (mean of repeats)", "CV-AUC signature min", "CV-AUC signature max",
             "CV-AUC age+sex", "CV-AUC signature+age+sex", "CV calibration slope (mean)", "CV calibration-in-the-large (mean)",
             "Fold-level AUC SD", "Median genes selected per fold"),
  value = c(mean(per_rep(oof)), min(per_rep(oof)), max(per_rep(oof)), mean(per_rep(oof_as)), mean(per_rep(oof_las)),
            mean(slope), mean(citl), sd(fold_auc), median(table(paste(sel$rep, sel$fold)))))
## pooled out-of-fold probability (average LP over repeats) -> Youden threshold (frozen)
lp_bar <- rowMeans(oof); p_bar <- plogis(lp_bar)
rc <- roc(y, p_bar, quiet = TRUE); ci_p <- ci.auc(rc)
thr <- coords(rc, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
cv_summary <- rbind(cv_summary, data.frame(metric = c("Pooled OOF AUC (averaged LP)", "Pooled OOF AUC lower 95% CI", "Pooled OOF AUC upper 95% CI",
  "Youden threshold (probability)", "OOF sensitivity at threshold", "OOF specificity at threshold"),
  value = c(ci_p[2], ci_p[1], ci_p[3], thr$threshold[1], thr$sensitivity[1], thr$specificity[1])))
freq <- as.data.frame(sort(table(sel$gene), decreasing = TRUE) / (R * K)); names(freq) <- c("gene", "selection_frequency")
dir.create("results", showWarnings = FALSE)
write.csv(cv_summary, "results/discovery_nestedcv_summary.csv", row.names = FALSE)
write.csv(data.frame(rep = 1:R, auc_signature = per_rep(oof), auc_age_sex = per_rep(oof_as), auc_sig_age_sex = per_rep(oof_las), cal_slope = slope), "results/discovery_cv_per_repeat.csv", row.names = FALSE)
write.csv(freq, "results/discovery_selection_frequency.csv", row.names = FALSE)
write.csv(data.frame(id = d$id, AD = y, oof_lp_mean = lp_bar, oof_as_mean = rowMeans(oof_as), oof_las_mean = rowMeans(oof_las)),
          "results/discovery_oof_predictions.csv", row.names = FALSE)
print(cv_summary, digits = 3)

write.csv(do.call(rbind, sel_log), "results/discovery_cv_fold_selection_log.csv", row.names = FALSE)
write.csv(data.frame(id = d$id, FOLDS), "results/discovery_cv_fold_assignment.csv", row.names = FALSE)
write.csv(data.frame(metric = c("Calibration slope, pooled OOF LP within repeat (mean of 10)", "Calibration slope, fold-wise (mean of 50 fold models on own held-out fold)", "Calibration slope, fold-wise median", "CITL, fold-wise mean", "Fold-level AUC mean", "Fold-level AUC SD"),
  value = c(mean(slope), mean(fold_slope), median(fold_slope), mean(fold_citl), mean(fold_auc), sd(fold_auc))), "results/discovery_cv_calibration_estimands.csv", row.names = FALSE)  # genes, lambda rule and penalty per outer fold
## --- final model on all discovery data (frozen) ---
set.seed(20261025)
final <- pipeline(seq_len(nrow(d)))
cat("Final genes (", length(final$genes), ", rule ", final$rule, ", penalty ", final$penalty, "):\n", paste(final$genes, collapse = ", "), "\n")
cf <- coef(final$fit)
frozen <- data.frame(term = c("Intercept", final$genes), coefficient = as.numeric(cf))
write.csv(frozen, "results/frozen_model_coefficients.csv", row.names = FALSE)
writeLines(c(paste0("threshold_probability,", thr$threshold[1]), paste0("penalty,", final$penalty), paste0("rule,", final$rule)),
           "results/frozen_model_settings.csv")
saveRDS(list(final = final, threshold = thr$threshold[1], genes = final$genes), "results/frozen_model.rds")
## --- elastic-net step of the final model: lambda path summary and coefficients (standardised predictors) ---
cvf <- final$cv; gl <- cvf$glmnet.fit; lam_used <- gl$lambda[max(which(gl$df <= 20))]
en <- data.frame(quantity = c("alpha", "candidate genes (step A)", "lambda.min", "non-zero genes at lambda.min", "lambda.1se", "non-zero genes at lambda.1se",
  "lambda used (least penalised with <= 20 genes)", "non-zero genes at lambda used", "CV deviance at lambda.min", "CV deviance at lambda used", "CV deviance SE at lambda used", "lambda path length"),
  value = c(0.5, length(final$candidates), cvf$lambda.min, cvf$nzero[cvf$lambda == cvf$lambda.min], cvf$lambda.1se, cvf$nzero[cvf$lambda == cvf$lambda.1se],
            lam_used, sum(coef(cvf, s = lam_used)[-1, 1] != 0), min(cvf$cvm), cvf$cvm[gl$lambda == lam_used], cvf$cvsd[gl$lambda == lam_used], length(gl$lambda)))
write.csv(en, "results/final_elastic_net_summary.csv", row.names = FALSE)
b_en <- coef(cvf, s = lam_used)[, 1]; b_en <- b_en[b_en != 0]
write.csv(data.frame(term = names(b_en), elastic_net_coefficient = unname(b_en)), "results/final_elastic_net_coefficients.csv", row.names = FALSE)
path <- data.frame(lambda = gl$lambda, nonzero = gl$df, cv_deviance = cvf$cvm, cv_se = cvf$cvsd); write.csv(path, "results/final_elastic_net_path.csv", row.names = FALSE)
## --- penalised refit: penalty trace (corrected AIC) on the original and an extended grid ---
dd0 <- data.frame(AD = y, X_all[, final$genes, drop = FALSE], check.names = FALSE); names(dd0) <- make.names(names(dd0))
f00 <- lrm(as.formula(paste("AD ~", paste(make.names(final$genes), collapse = " + "))), data = dd0, x = TRUE, y = TRUE)
ptr <- pentrace(f00, c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512))
write.csv(as.data.frame(ptr$results.all), "results/final_pentrace_extended_grid.csv", row.names = FALSE)
st <- final$fit$stats; write.csv(data.frame(statistic = names(st), value = unname(st)), "results/final_lrm_fit_statistics.csv", row.names = FALSE)
write.csv(data.frame(term = names(coef(f00)), unpenalised_ML = unname(coef(f00)), unpenalised_SE = unname(sqrt(diag(vcov(f00)))),
  penalised = unname(coef(final$fit)), penalised_SE = unname(sqrt(diag(vcov(final$fit))))), "results/final_lrm_unpenalised_vs_penalised.csv", row.names = FALSE)
cat("extended-grid optimum penalty:", ptr$penalty, "\n")

## --- full-pipeline bootstrap optimism (B = 500) ---
lp_full <- lp_new(final, X_all)
app_auc <- as.numeric(auc(y, lp_full, quiet = TRUE)); app_slope <- 1
B <- 500
library(parallel)
cl <- makeCluster(16); clusterSetRNGStream(cl, 20261025)
clusterExport(cl, c("d", "X_all", "y", "filter_genes", "select_en", "fit_lrm", "lp_new", "pipeline"))
clusterEvalQ(cl, {suppressPackageStartupMessages({library(glmnet); library(rms); library(pROC)}); NULL})
res <- parLapply(cl, 1:B, function(b) {
  bi <- sample(seq_len(nrow(d)), replace = TRUE)
  mb <- tryCatch(pipeline(bi), error = function(e) NULL); if (is.null(mb)) return(c(NA, NA, NA))
  lb <- lp_new(mb, X_all[bi, ]); lo <- lp_new(mb, X_all)
  c(as.numeric(auc(y[bi], lb, quiet = TRUE)) - as.numeric(auc(y, lo, quiet = TRUE)),
    1 - coef(glm(y ~ lo, binomial))[2], length(mb$genes))
})
stopCluster(cl)
res <- do.call(rbind, res); opt_auc <- res[, 1]; opt_slope <- res[, 2]
boot <- data.frame(metric = c("Apparent AUC", "Optimism (AUC)", "Optimism-corrected AUC", "Optimism-corrected calibration slope", "Bootstrap failures"),
                   value = c(app_auc, mean(opt_auc, na.rm = TRUE), app_auc - mean(opt_auc, na.rm = TRUE), 1 - mean(opt_slope, na.rm = TRUE), sum(is.na(opt_auc))))
write.csv(boot, "results/discovery_bootstrap_optimism.csv", row.names = FALSE)
print(boot, digits = 3)
