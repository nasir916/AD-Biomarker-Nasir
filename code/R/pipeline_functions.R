
## R/pipeline_functions.R — shared development pipeline (filter -> elastic net -> penalised lrm), parameterised
suppressPackageStartupMessages({library(glmnet); library(rms); library(pROC)})
META <- c("id","AD","dx","age","sex","batch","apoe")
read_cohort <- function(f) {
  x <- read.csv(gzfile(f), check.names = FALSE)
  if (!"AD" %in% names(x)) x$AD <- as.integer(x$dx %in% c("AD"))
  x$sex <- factor(tolower(x$sex), levels = c("female","male"))
  if ("batch" %in% names(x)) x$batch <- factor(x$batch)
  x }
filter_genes <- function(d, X, idx, top = 200, batch = TRUE) {
  f <- if (batch && "batch" %in% names(d)) ~ AD + age + sex + batch else ~ AD + age + sex
  D <- model.matrix(f, data = droplevels(d[idx, ])); Y <- X[idx, , drop = FALSE]
  q <- qr(D); B <- qr.coef(q, Y); R <- Y - D %*% B
  df <- nrow(D) - q$rank; s2 <- colSums(R^2) / df
  XtXi <- chol2inv(qr.R(q))[1:q$rank, 1:q$rank]; j <- which(colnames(D) == "AD")
  p <- 2 * pt(-abs(B[j, ] / sqrt(s2 * XtXi[j, j])), df)
  names(sort(p))[1:top] }
select_en <- function(X, y, idx, g, alpha = 0.5, cap = 20) {
  cv <- cv.glmnet(X[idx, g], y[idx], family = "binomial", alpha = alpha, nfolds = 10)
  nz <- function(l) { b <- coef(cv, s = l)[-1, 1]; names(b)[b != 0] }
  sel <- nz(cv$lambda.1se)
  if (is.finite(cap) && length(sel) > cap) { ok <- which(cv$glmnet.fit$df <= cap); sel <- nz(cv$glmnet.fit$lambda[max(ok)]) }
  if (length(sel) == 0) sel <- nz(cv$lambda.min)
  sel }
fit_lrm <- function(X, y, idx, g) {
  dd <- data.frame(AD = y[idx], X[idx, g, drop = FALSE], check.names = FALSE); names(dd) <- c("AD", make.names(g))
  f0 <- lrm(as.formula(paste("AD ~", paste(make.names(g), collapse = " + "))), data = dd, x = TRUE, y = TRUE)
  pen <- pentrace(f0, c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128))$penalty
  list(fit = if (pen > 0) update(f0, penalty = pen) else f0, penalty = pen, genes = g) }
lp_new <- function(m, newX) { nd <- as.data.frame(newX[, m$genes, drop = FALSE]); names(nd) <- make.names(m$genes); predict(m$fit, newdata = nd, type = "lp") }
pipeline <- function(d, X, idx, top = 200, alpha = 0.5, cap = 20, batch = TRUE) {
  g <- filter_genes(d, X, idx, top, batch); s <- select_en(X, d$AD, idx, g, alpha, cap); fit_lrm(X, d$AD, idx, s) }
strat_folds <- function(y, K, seed) { set.seed(seed); f <- integer(length(y)); for (cl in 0:1) { ii <- which(y == cl); f[ii] <- sample(rep(1:K, length.out = length(ii))) }; f }
auc_ci <- function(y, s) { r <- roc(y, s, quiet = TRUE, direction = "<"); as.numeric(ci.auc(r, method = "delong")) }  # lo, auc, hi
cal <- function(y, lp) c(intercept = unname(coef(glm(y ~ offset(lp), binomial))[1]), slope = unname(coef(glm(y ~ lp, binomial))[2]))
wilson <- function(x, n) { p <- x / n; z <- 1.96; c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n), (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }
