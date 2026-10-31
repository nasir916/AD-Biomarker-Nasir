
## 3_nomogram_tables.R — nomogram, points allocation and worked examples from the frozen model (discovery data only)
suppressPackageStartupMessages({library(rms)})
fm <- readRDS("results/frozen_model.rds"); fit <- fm$final$fit; G <- fm$genes
d <- read.csv(gzfile("data/GSE140829_z.csv.gz"), check.names = FALSE)
dd <- data.frame(AD = d$AD, d[, G], check.names = FALSE); names(dd) <- c("AD", make.names(G))
assign("ddist", datadist(dd), envir = .GlobalEnv); options(datadist = "ddist")
lims <- setNames(lapply(make.names(G), function(g) seq(-3, 3, by = 1)), make.names(G))
nom <- do.call(nomogram, c(list(fit, fun = plogis, funlabel = "Probability of AD",
        fun.at = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95), lp = FALSE, lp.at = seq(-5.5, 5, by = 0.5)), lims))
png("results/Fig3a_nomogram.png", width = 2400, height = 2600, res = 300)
par(mar = c(1, 1, 1, 1)); plot(nom, xfrac = 0.22, cex.axis = 0.62, cex.var = 0.8, lmgp = 0.2); dev.off()
## points per gene at z = -3..3
pts <- do.call(rbind, lapply(make.names(G), function(g) {
  x <- nom[[g]]; data.frame(gene = G[match(g, make.names(G))], z = x[[1]], points = round(x$points, 1)) }))
write.csv(pts, "results/Table4_points_allocation.csv", row.names = FALSE)
## total points -> probability
tp <- nom$`Probability of AD`; write.csv(data.frame(total_points = round(tp$x, 1), probability = tp$x.real),
  "results/Table4b_points_to_probability.csv", row.names = FALSE)
## worked examples: the discovery AD case and control closest to the median predicted risk of their group
lp <- predict(fit, dd, type = "lp"); pick <- function(cl) { i <- which(d$AD == cl); i[which.min(abs(lp[i] - median(lp[i])))] }
ex <- c(pick(1), pick(0))
pf <- function(i) { s <- sapply(make.names(G), function(g) approx(nom[[g]][[1]], nom[[g]]$points, xout = dd[i, g], rule = 2)$y); s }
tab <- data.frame(gene = G, example1_z = round(unlist(dd[ex[1], make.names(G)]), 2), example1_points = round(pf(ex[1]), 1),
                  example2_z = round(unlist(dd[ex[2], make.names(G)]), 2), example2_points = round(pf(ex[2]), 1))
tot <- colSums(tab[, c("example1_points", "example2_points")])
tab <- rbind(tab, data.frame(gene = "Total points", example1_z = NA, example1_points = round(tot[1], 1), example2_z = NA, example2_points = round(tot[2], 1)),
             data.frame(gene = "Predicted probability (model)", example1_z = NA, example1_points = round(plogis(lp[ex[1]]), 3), example2_z = NA, example2_points = round(plogis(lp[ex[2]]), 3)),
             data.frame(gene = "Observed diagnosis", example1_z = NA, example1_points = d$AD[ex[1]], example2_z = NA, example2_points = d$AD[ex[2]]))
write.csv(tab, "results/Table5_worked_examples.csv", row.names = FALSE)
cat("max attainable total points:", sum(sapply(make.names(G), function(g) max(nom[[g]]$points))), " probability axis max points:", max(tp$x), " examples:", ex, " totals:", tot, " probs:", plogis(lp[ex]), "\n")
## coefficient table with OR per SD and SE
se <- sqrt(diag(vcov(fit)))
ct <- data.frame(term = names(coef(fit)), beta = coef(fit), SE = se, OR_per_SD = exp(coef(fit)),
                 OR_lo = exp(coef(fit) - 1.96 * se), OR_hi = exp(coef(fit) + 1.96 * se))
ct$term[-1] <- G; write.csv(ct, "results/Table3_coefficients.csv", row.names = FALSE)
