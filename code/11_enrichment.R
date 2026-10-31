
## 11_enrichment.R — genome-wide ranked gene-set enrichment of AD-associated expression against blood transcription modules
## (Li et al. 2014, distributed in the tmod R package; exported once to data/genesets/Li2014_BTM.gmt and read from there).
## Primary test: competitive gene-set test on t-statistics with a fixed inter-gene correlation (cameraPR algorithm of limma,
## re-implemented below; limma is not installable on this platform), with a rank-based (correlation-adjusted Wilcoxon) check.
source("R/pipeline_functions.R"); OUT <- "results/sensitivity/"; dir.create("data/genesets", showWarnings = FALSE)
gmt <- "data/genesets/Li2014_BTM.gmt"
if (!file.exists(gmt)) { e1 <- new.env(); lib <- c(normalizePath("../Rlib", mustWork = FALSE), normalizePath("Rlib", mustWork = FALSE))
  f <- file.path(lib, "tmod", "data", "tmod.rda"); f <- f[file.exists(f)][1]; load(f, envir = e1); tm <- e1$tmod
  li <- tm$gs[tm$gs$SourceID == "LI", ]; idx <- match(li$ID, tm$gs$ID)
  lines <- sapply(seq_len(nrow(li)), function(i) paste(c(li$ID[i], li$Title[i], tm$gv[tm$gs2gv[[idx[i]]]]), collapse = "\t")); writeLines(lines, gmt)
  writeLines(c(paste("source: tmod R package, version", as.character(read.dcf(file.path(dirname(dirname(f)), "DESCRIPTION"))[, "Version"])), "modules: Li et al. 2014 (Nat Immunol 15:195-204), SourceID LI"), "data/genesets/Li2014_BTM_source.txt") }
GS <- lapply(strsplit(readLines(gmt), "\t"), function(x) list(id = x[1], title = x[2], genes = unique(x[-(1:2)])))
## t-statistics
tstat <- function(df, X, f) { D <- model.matrix(f, data = droplevels(df)); q <- qr(D); B <- qr.coef(q, X); R <- X - D %*% B; s2 <- colSums(R^2) / (nrow(D) - q$rank)
  XtXi <- chol2inv(qr.R(q)); j <- which(colnames(D) == "AD"); B[j, ] / sqrt(s2 * XtXi[j, j]) }
d <- read_cohort("data/GSE140829_z.csv.gz"); genes <- setdiff(names(d), META); Xd <- as.matrix(d[, genes])
TT <- data.frame(gene = genes, disc_batch_adj = tstat(d, Xd, ~ AD + age + sex + batch), disc_age_sex = tstat(d, Xd, ~ AD + age + sex))
for (cn in c("GSE63061", "GSE63060")) { w <- read_cohort(sprintf("data/%s_retained_fixedref.csv.gz", cn)); w$AD <- as.integer(w$dx == "AD")
  inf <- read.csv(gzfile(sprintf("data/%s_inclSexDiscordant_z.csv.gz", cn)), check.names = FALSE)[, c("id", "sex_inferred")]; w$sex <- factor(inf$sex_inferred[match(w$id, inf$id)])
  TT[[cn]] <- tstat(w, as.matrix(w[, genes]), ~ AD + age + sex) }
write.csv(TT, paste0(OUT, "AK_tstats_for_enrichment.csv"), row.names = FALSE)
cameraPR1 <- function(stat, idx, cor = 0.01) { G <- length(stat); m <- length(idx); ms <- mean(stat); vs <- var(stat)
  delta <- G / (G - m) * (mean(stat[idx]) - ms); vp <- ((G - 1) * vs - delta^2 * m * (G - m) / G) / (G - 2); vif <- 1 + (m - 1) * cor
  tt <- delta / sqrt(vp * (vif / m + 1 / (G - m))); c(t = tt, p = 2 * pt(-abs(tt), df = G - 2)) }
ranksum_cor <- function(stat, idx, cor = 0.01) { n <- length(stat); r <- rank(stat); r1 <- r[idx]; n1 <- length(r1); n2 <- n - n1
  U <- n1 * n2 + n1 * (n1 + 1) / 2 - sum(r1); mu <- n1 * n2 / 2
  s2 <- (asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 1) + asin(cor / 2) * n1 * (n1 - 1) * n2 * (n2 - 1) + asin((cor + 1) / 2) * n1 * (n1 - 1) * n2) / 2 / pi
  z <- (U - mu) / sqrt(s2); 2 * pnorm(-abs(z)) }
cls <- function(t) { t <- tolower(t)
  if (grepl("monocyt|neutrophil|myeloid|inflammat|tlr|innate|complement|dendritic", t)) "myeloid/inflammatory"
  else if (grepl("b cell|t cell|nk cell|natural killer|lymph|cytotoxic|plasma cell|immunoglobulin", t)) "lymphoid" else "other" }
RES <- list()
for (col in c("disc_batch_adj", "disc_age_sex", "GSE63061", "GSE63060")) { st <- setNames(TT[[col]], TT$gene)
  R <- do.call(rbind, lapply(GS, function(s) { idx <- which(names(st) %in% s$genes); if (length(idx) < 5 || length(idx) > 500) return(NULL)
    cp <- cameraPR1(st, idx); data.frame(module = s$id, title = s$title, class = cls(s$title), n_genes = length(idx), t = cp["t"], direction = ifelse(cp["t"] > 0, "up", "down"),
      p = cp["p"], p_rank = ranksum_cor(st, idx)) }))
  R$FDR <- p.adjust(R$p, "BH"); R$FDR_rank <- p.adjust(R$p_rank, "BH"); names(R)[5:10] <- paste0(col, "_", names(R)[5:10]); RES[[col]] <- R }
E <- Reduce(function(a, b) merge(a, b, by = c("module", "title", "class", "n_genes")), RES)
E <- E[order(E$disc_batch_adj_p), ]; write.csv(E, paste0(OUT, "AK_enrichment_Li_BTM.csv"), row.names = FALSE)
## pre-specified summary: modules with FDR < 0.05 in discovery (batch-adjusted) and same direction with p < 0.05 in each validation cohort
sig <- E$disc_batch_adj_FDR < 0.05
S <- do.call(rbind, lapply(c("myeloid/inflammatory", "lymphoid", "other"), function(k) { i <- E$class == k
  data.frame(class = k, modules_tested = sum(i), disc_FDR05_up = sum(sig & i & E$disc_batch_adj_direction == "up"), disc_FDR05_down = sum(sig & i & E$disc_batch_adj_direction == "down"),
    replicated_both_validation = sum(sig & i & E$GSE63061_direction == E$disc_batch_adj_direction & E$GSE63060_direction == E$disc_batch_adj_direction & E$GSE63061_p < 0.05 & E$GSE63060_p < 0.05),
    up_FDR05_in_at_least_2_of_3 = sum(i & ((E$disc_batch_adj_FDR < 0.05 & E$disc_batch_adj_direction == "up") + (E$GSE63061_FDR < 0.05 & E$GSE63061_direction == "up") + (E$GSE63060_FDR < 0.05 & E$GSE63060_direction == "up")) >= 2),
    down_FDR05_in_at_least_2_of_3 = sum(i & ((E$disc_batch_adj_FDR < 0.05 & E$disc_batch_adj_direction == "down") + (E$GSE63061_FDR < 0.05 & E$GSE63061_direction == "down") + (E$GSE63060_FDR < 0.05 & E$GSE63060_direction == "down")) >= 2)) }))
write.csv(S, paste0(OUT, "AK_enrichment_summary.csv"), row.names = FALSE); print(S); print(head(E[, c(1:4, 5, 7, 9, 11, 13, 17, 19)], 25), digits = 3)
