
"""4_figures.py — Figures 1-5 and Supplementary Figures S1 (nomogram) and S2 from results/ CSVs."""
import os, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
os.makedirs("figs", exist_ok=True); R_ = "results/"; RV = R_ + "sensitivity/"
mpl.rcParams.update({"font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
                     "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
C = {"GSE63061": "#1f5fa8", "GSE63060": "#d07c1d", "disc": "#444444"}
def roc_curve(y, s):
    o = np.argsort(-np.asarray(s)); y = np.asarray(y)[o]; tps = np.cumsum(y); fps = np.cumsum(1 - y)
    return np.r_[0, fps / fps[-1]], np.r_[0, tps / tps[-1]]
def letter(ax, l): ax.text(-0.14, 1.06, l, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
PA = pd.read_csv(RV + "AA_primary_external_performance.csv").set_index("cohort")
# ---------- Fig 1 flowchart ----------
fig, ax = plt.subplots(figsize=(7.4, 5.6)); ax.axis("off"); ax.set_xlim(-0.2, 10.2); ax.set_ylim(0, 12.6)
def box(x, y, w, h, t, fc="#f2f2f2", fs=6.6):
    ax.add_patch(mpl.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec="#555", lw=0.8))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs)
def arr(x1, y1, x2, y2): ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#555"))
box(1.0, 11.3, 8.0, 1.0, "GEO search: 197 human AD blood-related series screened\n3 Illumina HumanHT-12 whole-blood series with ≥ 40 AD and ≥ 40 CN")
box(0.05, 8.2, 3.15, 2.5, "Discovery: GSE140829\nHumanHT-12 v4, USA\nAD 204 / CN 249 before QC\nQC: 2 sex-discordant, 9 outliers\nAnalysed: AD 200 / CN 242", fc="#e8e8e8")
box(3.42, 8.2, 3.15, 2.5, "Validation: GSE63061\nAddNeuroMed batch 2, HumanHT-12 v4\nAD 139 / CN 134 before QC\nQC: 5 outliers removed\nAnalysed: AD 136 / CN 132\n(4 sex-discordant retained)", fc="#dfe9f5")
box(6.8, 8.2, 3.15, 2.5, "Validation: GSE63060\nAddNeuroMed batch 1, HumanHT-12 v3\nAD 145 / CN 104 before QC\nQC: 2 outliers removed\nAnalysed: AD 143 / CN 104\n(49 sex-discordant retained)", fc="#f8ead8")
for x in (1.62, 5.0, 8.37): arr(5, 11.3, x, 10.75)
box(0.05, 5.0, 3.15, 2.6, "Shared universe: 9,804 genes\nPipeline inside 5×10 nested CV:\n(A) adjusted filter, top 200\n(B) elastic net, ≤ 20 genes\n(C) penalised logistic regression", fc="#e8e8e8")
arr(1.62, 8.2, 1.62, 7.65)
box(0.05, 2.9, 3.15, 1.6, "Final 17-gene model\n(full discovery cohort)\nthreshold p = 0.451", fc="#d9d9d9"); arr(1.62, 5.0, 1.62, 4.55)
box(3.7, 4.2, 6.2, 3.4, "Model applied without modification\nfixed-reference scaling (primary); within-cohort z (sensitivity)\n\nAUC, calibration, Brier score; added value over age + sex\nrandom-signature benchmarks (4 types)\nbatch, cell composition, age, sex, probe and site analyses\nsex-discordant samples excluded (sensitivity)", fc="#ffffff")
arr(3.25, 3.7, 3.65, 4.6); arr(5.0, 8.2, 6.0, 7.65); arr(8.37, 8.2, 7.6, 7.65)
box(3.7, 0.2, 6.2, 3.4, "Specificity (cross-fitted fold models)\nGSE140830: 542 arrays → QC (3 sex-discordant, 14 outliers) → 526\n240 control arrays duplicating GSE140829 arrays removed\nAnalysed: bvFTD 80, PSP 51, nfvPPA 46, svPPA 41, CBS 35,\n33 independent controls\nGSE140829 MCI: 134 → QC (7 outliers) → 127", fc="#eef5ea")
box(0.05, 0.6, 3.15, 1.8, "Gene-set enrichment\n(blood transcription modules)\ndiscovery and both validation cohorts", fc="#f2f2f2"); arr(1.62, 2.9, 1.62, 2.45)
fig.savefig("figs/Fig1_flowchart.png", dpi=300, bbox_inches="tight"); plt.close(fig)
# ---------- Fig 2 discovery ----------
rep = pd.read_csv(R_ + "discovery_cv_per_repeat.csv"); cvs = pd.read_csv(R_ + "discovery_nestedcv_summary.csv").set_index("metric").value
fz = pd.read_csv(R_ + "frozen_model_coefficients.csv"); fr = pd.read_csv(R_ + "discovery_selection_frequency.csv")
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw=dict(width_ratios=[1, 1.4], wspace=0.45))
ax = axs[0]; rng = np.random.default_rng(1)
for i, (col, lab) in enumerate([("auc_age_sex", "Age + sex"), ("auc_signature", "17-gene\nsignature"), ("auc_sig_age_sex", "Signature +\nage + sex")]):
    v = rep[col].values; ax.scatter(i + rng.uniform(-.08, .08, len(v)), v, s=10, color=C["disc"], alpha=.6, zorder=3)
    ax.hlines(v.mean(), i - 0.28, i + 0.28, color="k", lw=1.6); ax.text(i + 0.32, v.mean(), f"{v.mean():.2f}", ha="left", va="center", fontsize=6)
ax.axhline(0.5, ls=":", color="#999", lw=0.7); ax.set_xticks(range(3)); ax.set_xticklabels(["Age + sex", "17-gene\nsignature", "Signature +\nage + sex"])
ax.set_ylabel("Nested-CV AUC, GSE140829 (10 repeats)"); ax.set_ylim(0.4, 0.72); ax.set_xlim(-0.5, 2.9); ax.set_title("Internal validation", loc="left"); letter(ax, "a")
ax = axs[1]; g = fz.iloc[1:].copy(); g["freq"] = g.term.map(dict(zip(fr.gene, fr.selection_frequency))); g = g.sort_values("coefficient").reset_index(drop=True)
ax.barh(g.term, g.coefficient, color=np.where(g.coefficient > 0, "#b2182b", "#2166ac"), height=0.65)
for yy, f in enumerate(g.freq): ax.text(0.262, yy, f"{f:.0%}", va="center", ha="right", fontsize=6, color="#555")
ax.text(0.262, len(g) - 0.4, "% of CV folds", fontsize=6, ha="right", va="bottom", color="#555")
ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("Penalised log-odds of AD per SD of expression"); ax.set_xlim(-0.14, 0.27)
for t in ax.get_yticklabels(): t.set_style("italic")
ax.set_title("Final 17-gene model", loc="left"); letter(ax, "b")
fig.savefig("figs/Fig2_discovery.png", dpi=300, bbox_inches="tight"); plt.close(fig)
# ---------- Fig 3 external validation (primary analysis) ----------
pr = pd.read_csv(RV + "AA_primary_external_predictions.csv"); cur = pd.read_csv(RV + "AA_primary_loess_curves.csv")
BB = pd.read_csv(RV + "AB_primary_random_benchmarks_draws.csv"); TT = pd.read_csv(RV + "AK_tstats_for_enrichment.csv").set_index("gene"); G = fz.term[1:].tolist()
EN = pd.read_csv(RV + "AK_enrichment_Li_BTM.csv")
fig = plt.figure(figsize=(7.8, 5.8)); gs = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.5, width_ratios=[1, 1, 1.0])
ax = fig.add_subplot(gs[0, 0])
for c in ["GSE63061", "GSE63060"]:
    s = pr[pr.cohort == c]; fpr, tpr = roc_curve(s.AD, s.lp)
    ax.plot(fpr, tpr, color=C[c], lw=1.3, label=f"{c}: {PA.loc[c,'AUC']:.2f} ({PA.loc[c,'AUC_lo']:.2f}–{PA.loc[c,'AUC_hi']:.2f})")
ax.plot([0, 1], [0, 1], ls=":", color="#999", lw=0.7); ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
ax.legend(loc="lower right", frameon=False, fontsize=5.5, title="AUC (95% CI)", title_fontsize=5.5); ax.set_title("Discrimination", loc="left"); letter(ax, "a")
ax = fig.add_subplot(gs[0, 1])
for c in ["GSE63061", "GSE63060"]:
    s = pr[pr.cohort == c].copy(); s["bin"] = pd.qcut(s.p, 5, labels=False, duplicates="drop"); b = s.groupby("bin").agg(p=("p", "mean"), o=("AD", "mean"), n=("AD", "size"))
    se = np.sqrt(b.o * (1 - b.o) / b.n); ax.errorbar(b.p, b.o, yerr=1.96 * se, fmt="o", ms=3, color=C[c], lw=0.8, capsize=1.5, label=c)
    cc = cur[cur.cohort == c]; ax.plot(cc.predicted, cc.observed_loess, color=C[c], lw=1.1)
ax.plot([0, 1], [0, 1], ls=":", color="#999", lw=0.7); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed proportion AD"); ax.legend(loc="upper left", frameon=False, fontsize=5.5)
ax.set_title("Calibration", loc="left"); letter(ax, "b")
ax = fig.add_subplot(gs[0, 2])
order = ["expression-matched", "composition-matched", "step-A candidate pool", "strong-association pool"]
labs = ["Expr.-\nmatch", "Comp.-\nmatch", "Step-A\npool", "Strong\npool"]
for j, c in enumerate(["GSE63061", "GSE63060"]):
    col = "a61" if c == "GSE63061" else "a60"; pos = np.array([0, 1, 2.6, 3.6]) + (j - 0.5) * 0.36
    bp = ax.boxplot([BB.loc[BB.type == t, col].values for t in order], positions=pos, widths=0.3, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=0.8))
    for patch in bp["boxes"]: patch.set(facecolor=C[c], alpha=0.45, lw=0.6)
    ax.scatter(pos, [PA.loc[c, "AUC"]] * 4, marker="D", s=12, color=C[c], zorder=4, label=c)
ax.axvline(1.8, color="#bbb", lw=0.6); ax.text(0.5, 0.81, "whole universe", ha="center", fontsize=5.3, color="#666"); ax.text(3.1, 0.81, "discovery pools", ha="center", fontsize=5.3, color="#666")
ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(labs, fontsize=5.0); ax.set_ylim(0.4, 0.83); ax.set_ylabel("Validation AUC")
ax.legend(frameon=False, fontsize=5, loc="lower right"); ax.set_title("Random 17-gene sets (◆ model)", loc="left", pad=9); letter(ax, "c")
for k, c in enumerate(["GSE63061", "GSE63060"]):
    ax = fig.add_subplot(gs[1, k]); ax.scatter(TT["disc_age_sex"], TT[c], s=2, color="#cccccc", rasterized=True, label="All 9,804 genes")
    ax.scatter(TT.loc[G, "disc_age_sex"], TT.loc[G, c], s=12, color=C[c], zorder=3, label="Signature genes")
    rho = TT["disc_age_sex"].corr(TT[c], method="spearman"); ax.axhline(0, color="#999", lw=0.5); ax.axvline(0, color="#999", lw=0.5)
    ax.set_xlabel("AD t-statistic, GSE140829"); ax.set_ylabel(f"AD t-statistic, {c}"); ax.text(0.02, 0.98, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6)
    ax.legend(frameon=False, fontsize=5.3, loc="lower right"); ax.set_title(f"Genome-wide AD effects, {c}", loc="left"); letter(ax, "de"[k])
ax = fig.add_subplot(gs[1, 2])
top = EN[EN["class"] != "other"].nsmallest(10, "disc_batch_adj_p").iloc[::-1]
for j, (col, mk_, cc) in enumerate([("disc_batch_adj_t", "s", "#444444"), ("GSE63061_t", "o", C["GSE63061"]), ("GSE63060_t", "^", C["GSE63060"])]):
    ax.scatter(top[col], np.arange(len(top)) + (j - 1) * 0.22, marker=mk_, s=12, color=cc, label=["GSE140829", "GSE63061", "GSE63060"][j])
ax.axvline(0, color="#999", lw=0.6); ax.set_yticks(np.arange(len(top))); ax.yaxis.tick_right(); ax.spines["right"].set_visible(True); ax.spines["left"].set_visible(False); ax.set_yticklabels([t[:30] for t in top.title], fontsize=5)
ax.set_xlabel("Module statistic (camera t)"); ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3); ax.set_title("Blood transcription modules", loc="left"); letter(ax, "f")
fig.savefig("figs/Fig3_validation.png", dpi=300, bbox_inches="tight"); plt.close(fig)
# ---------- Fig 4 drivers ----------
gr = pd.read_csv(RV + "A_robustness_grid.csv"); sb = pd.read_csv(RV + "R_development_subsampling_draws.csv")
lb = pd.read_csv(RV + "C_lobo_summary.csv").set_index("metric").value; lb2 = pd.read_csv(RV + "C_lobo_batchstandardised_summary.csv").set_index("metric").value
Qb = pd.read_csv(RV + "Q_batch_residualised_development.csv").set_index("analysis"); AI = pd.read_csv(RV + "AI_lobo_heldout_self_standardised.csv").set_index("analysis")
GC = pd.read_csv(RV + "AG_primary_composition.csv").set_index("cohort")
lab = {"base": "Main configuration", "alpha0.25": "α = 0.25", "alpha0.75": "α = 0.75", "lasso": "Lasso (α = 1)", "cap10": "≤ 10 genes", "cap15": "≤ 15 genes",
       "cap30": "≤ 30 genes", "nocap": "No gene cap", "top100": "Filter top 100", "top500": "Filter top 500", "ctlref": "Control-referenced scaling", "fixedref": "Fixed-reference pipeline"}
fig, axs = plt.subplots(1, 3, figsize=(7.8, 3.9), gridspec_kw=dict(width_ratios=[1.25, 1.05, 0.95], wspace=1.1))
ax = axs[0]; yy = np.arange(len(gr))[::-1]
for c, off in [("61", 0.15), ("60", -0.15)]:
    col = C["GSE63061"] if c == "61" else C["GSE63060"]
    ax.errorbar(gr[f"AUC{c}"], yy + off, xerr=[gr[f"AUC{c}"] - gr[f"AUC{c}_lo"], gr[f"AUC{c}_hi"] - gr[f"AUC{c}"]], fmt="o", ms=3, color=col, lw=0.8, label="GSE63061" if c == "61" else "GSE63060")
q = sb.AUC61.quantile([0.025, 0.5, 0.975]).values
ax.axvspan(q[0], q[2], color=C["GSE63061"], alpha=0.10, lw=0); ax.axvline(q[1], color=C["GSE63061"], lw=0.7, ls="--")
ax.set_yticks(yy); ax.set_yticklabels([lab[c] for c in gr.config]); ax.axvline(0.5, color="#999", lw=0.6, ls=":")
ax.set_xlabel("Validation AUC (95% CI)\nwithin-cohort scaling"); ax.set_xlim(0.45, 0.85); ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=5.5)
ax.set_title("Analytic choices", loc="left"); letter(ax, "a")
ax = axs[1]
rows = [("Random folds", cvs["CV-AUC signature (mean of repeats)"], cvs["CV-AUC signature min"], cvs["CV-AUC signature max"], "s", "#444444", "range"),
        ("Random folds,\nbatch-residualised", Qb.loc["Random-fold CV, mean per-repeat AUC", "estimate"], Qb.loc["Random-fold CV, mean per-repeat AUC", "lo"], Qb.loc["Random-fold CV, mean per-repeat AUC", "hi"], "s", "#7b3294", "range"),
        ("LOBO", lb["LOBO pooled AUC"], lb["LOBO pooled AUC lo"], lb["LOBO pooled AUC hi"], "o", "#444444", "ci"),
        ("LOBO, held-out batch\nstandardised on itself", AI.iloc[0]["pooled_AUC"], AI.iloc[0]["lo"], AI.iloc[0]["hi"], "o", "#2c7fb8", "ci"),
        ("LOBO, all batches\nstandardised", lb2["LOBO (within-batch standardised) pooled AUC"], lb2["lo"], lb2["hi"], "o", "#444444", "ci"),
        ("LOBO,\nbatch-residualised", Qb.loc["LOBO pooled AUC", "estimate"], Qb.loc["LOBO pooled AUC", "lo"], Qb.loc["LOBO pooled AUC", "hi"], "o", "#7b3294", "ci")]
for i, (n, v, lo, hi, mk_, col, kind) in enumerate(rows):
    ax.errorbar(v, len(rows) - 1 - i, xerr=[[v - lo], [hi - v]], fmt=mk_, ms=3.5, color=col, lw=0.9, ls="-" if kind == "ci" else "--", capsize=2 if kind == "range" else 0)
ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6); ax.axvline(0.5, color="#999", lw=0.6, ls=":")
ax.set_xlim(0.45, 0.75); ax.set_xlabel("Discovery AUC\n■ range over 10 repeats   ● DeLong 95% CI"); ax.set_title("Processing batch (GSE140829)", loc="left"); letter(ax, "b")
ax = axs[2]
cr = [("17-gene signature", "sig"), ("5 cell scores", "cells"), ("Cell scores +\nsignature", "cellsig"), ("Age, sex, 5 cell\nscores", "comp"), ("Age, sex, cell\nscores + signature", "compsig")]
for j, c in enumerate(["GSE63061", "GSE63060"]):
    v = {"sig": (PA.loc[c, "AUC"], PA.loc[c, "AUC_lo"], PA.loc[c, "AUC_hi"]), "cells": (GC.loc[c, "AUC_cells_model_from_discovery"], GC.loc[c, "AUC_cells_lo"], GC.loc[c, "AUC_cells_hi"]),
         "cellsig": (GC.loc[c, "AUC_cells_plus_sig_from_discovery"], np.nan, np.nan), "comp": (GC.loc[c, "AUC_comp_model_from_discovery"], np.nan, np.nan), "compsig": (GC.loc[c, "AUC_comp_plus_sig_from_discovery"], np.nan, np.nan)}
    for i, (n, k) in enumerate(cr):
        a, lo, hi = v[k]; yv = len(cr) - 1 - i + (0.15 if j == 0 else -0.15)
        ax.errorbar(a, yv, xerr=None if np.isnan(lo) else [[a - lo], [hi - a]], fmt="o", ms=3, color=C[c], lw=0.8, label=c if i == 0 else None)
ax.set_yticks(range(len(cr))[::-1]); ax.set_yticklabels([r[0] for r in cr], fontsize=5.6); ax.axvline(0.5, color="#999", lw=0.6, ls=":")
ax.set_xlim(0.45, 0.8); ax.set_xlabel("Validation AUC\n(models developed in discovery)"); ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax.set_title("Cell composition", loc="left"); letter(ax, "c")
fig.savefig("figs/Fig4_drivers.png", dpi=300, bbox_inches="tight"); plt.close(fig)
# ---------- Fig 5 specificity ----------
SA = pd.read_csv(RV + "AH_specificity_age_sex_adjusted.csv").set_index("group"); GT = pd.read_csv(RV + "AH_specificity_group_age_sex.csv").set_index("group")
sp = pd.read_csv(RV + "K_specificity_matched_ensemble.csv"); sp = sp[sp.reference == "Control (GSE140829)"].set_index("case")
order = ["AD (GSE140829)", "PSP", "nfvPPA", "CBS", "bvFTD", "svPPA", "MCI (GSE140829)", "Independent control (GSE140830)"]
nm = lambda g: g.replace(" (GSE140829)", "").replace("Independent control (GSE140830)", "Independent controls")
ypos = np.array([8, 6.6, 5.8, 5.0, 4.2, 3.4, 2.0, 0.8]); cols = ["#b2182b"] + ["#555555"] * 5 + ["#1a9850", "#999999"]
fig, axs = plt.subplots(1, 3, figsize=(7.8, 3.5), gridspec_kw=dict(width_ratios=[1.15, 0.9, 0.9], wspace=0.25))
ax = axs[0]
for g, yv, col in zip(order, ypos, cols):
    r = sp.loc[g]; ax.errorbar(r.AUC, yv + 0.15, xerr=[[r.AUC - r.lo], [r.hi - r.AUC]], fmt="o", ms=3.2, color=col, lw=0.8)
    ax.scatter(SA.loc[g, "AUC_age_matched_pairs"], yv - 0.17, marker="o", facecolor="white", edgecolor=col, s=14, zorder=3)
ax.set_yticks(ypos); ax.set_yticklabels([f"{nm(g)} (n={int(GT.loc[g,'n'])}; {GT.loc[g,'age_mean']:.0f} y)" for g in order], fontsize=5.8)
ax.axvline(0.5, color="#999", lw=0.6, ls=":"); ax.set_xlim(0.35, 0.85); ax.set_ylim(0.2, 8.6)
ax.set_xlabel("AUC vs discovery controls (73.7 y)\n● unadjusted (95% CI)   ○ age-matched pairs"); ax.set_title("Score by group", loc="left"); letter(ax, "a")
ax = axs[1]
for g, yv, col in zip(order, ypos, cols):
    r = SA.loc[g]; ax.errorbar(r.OR_per_SD_age_sex_adjusted, yv, xerr=[[r.OR_per_SD_age_sex_adjusted - r.OR_lo], [r.OR_hi - r.OR_per_SD_age_sex_adjusted]], fmt="s", ms=3, color=col, lw=0.8)
ax.set_xscale("log"); ax.axvline(1, color="#999", lw=0.6, ls=":"); ax.set_yticks(ypos); ax.set_yticklabels([]); ax.set_ylim(0.2, 8.6)
ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(["0.5", "1", "2", "4"]); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_xlabel("OR per SD of score vs discovery controls\nadjusted for age and sex (95% CI)")
ax.set_title("Age- and sex-adjusted", loc="left"); letter(ax, "b")
ax = axs[2]
for g, yv, col in zip(order[1:], ypos[1:], cols[1:]):
    r = SA.loc[g]; ax.errorbar(r.dAUC_residualised, yv, xerr=[[r.dAUC_residualised - r.dAUC_lo], [r.dAUC_hi - r.dAUC_residualised]], fmt="o", ms=3, color=col, lw=0.8)
ax.axvline(0, color="#999", lw=0.6, ls=":"); ax.set_yticks(ypos); ax.set_yticklabels([]); ax.set_ylim(0.2, 8.6)
ax.set_xlabel("ΔAUC: (AD vs CN) − (group vs CN)\nage/sex-residualised score, bootstrap 95% CI"); ax.set_title("Contrast with AD", loc="left"); letter(ax, "c")
fig.savefig("figs/Fig5_specificity.png", dpi=300, bbox_inches="tight"); plt.close(fig)
# ---------- Fig S1 nomogram (drawn from the points table) ----------
pts = pd.read_csv(R_ + "Table4_points_allocation.csv"); p2p = pd.read_csv(R_ + "Table4b_points_to_probability.csv")
genes = fz.term[1:].tolist(); maxpts = pts.groupby("gene").points.max().sum()
a1, a0 = np.polyfit(p2p.total_points, np.log(p2p.probability / (1 - p2p.probability)), 1)   # logit(p) = a0 + a1 * total points
fig, ax = plt.subplots(figsize=(7.2, 7.6)); ax.axis("off"); W = 100.0
nrow = len(genes) + 3; ytop = nrow
def axisline(yv, ticks, labels, label, scale=1.0, fs=5.8):
    ax.plot([min(ticks) * scale, max(ticks) * scale], [yv, yv], color="k", lw=0.7)
    for t, l in zip(ticks, labels): ax.plot([t * scale, t * scale], [yv, yv + 0.15], color="k", lw=0.6); ax.text(t * scale, yv - 0.12, l, ha="center", va="top", fontsize=fs)
    ax.text(-4, yv, label, ha="right", va="center", fontsize=6.5)
axisline(ytop, np.arange(0, 101, 10), [str(v) for v in range(0, 101, 10)], "Points")
for i, g in enumerate(genes):
    s = pts[pts.gene == g].set_index("z").points; yv = ytop - 1 - i
    zl = [-3, 0, 3] if s.max() < 25 else list(range(-3, 4))
    ax.plot([s.min(), s.max()], [yv, yv], color="k", lw=0.7)
    for z in zl: ax.plot([s[z], s[z]], [yv, yv + 0.15], color="k", lw=0.6); ax.text(s[z], yv - 0.12, str(z), ha="center", va="top", fontsize=5.5)
    ax.text(-4, yv, g, ha="right", va="center", fontsize=6.5, style="italic")
yT = ytop - 1 - len(genes); sc = W / maxpts; tk = np.arange(0, 701, 100)
axisline(yT, list(tk) + [maxpts], [str(int(v)) for v in tk] + [f"{maxpts:.0f}"], "Total points", scale=sc)
pr_t = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]; Tp = [(np.log(p / (1 - p)) - a0) / a1 for p in pr_t]
keep = [(p, t) for p, t in zip(pr_t, Tp) if 0 <= t <= maxpts]
yP = yT - 1.1; ax.plot([keep[0][1] * sc, keep[-1][1] * sc], [yP, yP], color="k", lw=0.7)
for p, t in keep: ax.plot([t * sc, t * sc], [yP, yP + 0.15], color="k", lw=0.6); ax.text(t * sc, yP - 0.12, f"{p:g}", ha="center", va="top", fontsize=5.5)
ax.text(-4, yP, "Probability of AD", ha="right", va="center", fontsize=6.5)
ax.text(50, ytop + 1.0, "NOT FOR CLINICAL USE — documents the model only. Inputs are fixed-reference z-scores; probabilities reflect the discovery prevalence (45%) and require local recalibration.",
        ha="center", va="bottom", fontsize=5.8, color="#b2182b", wrap=True)
ax.set_xlim(-30, 104); ax.set_ylim(yP - 1, ytop + 1.8)
fig.savefig("figs/FigS1_nomogram.png", dpi=300, bbox_inches="tight"); plt.close(fig)
print("figures written")
