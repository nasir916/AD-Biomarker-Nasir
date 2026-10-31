
"""9a_descriptive_tables.py — Table 1 (analysed cohorts, primary analysis sets, and specificity groups), Table 2 (signature genes),
agreement of genome-wide AD effects between cohorts with a direction-concordance null, checksums and package versions.
Reads data/ and results/ only."""
import json, hashlib, glob, os, sys, numpy as np, pandas as pd, scipy, matplotlib, statsmodels
from scipy import stats
R_ = "results/"; RV = R_ + "sensitivity/"
eth = pd.read_csv(RV + "P_ethnicity_AddNeuroMed.csv").fillna(0).set_index("cohort")
def white(c):
    r = eth.loc[c]; tot = r.sum(); nonw = sum(r.get(k, 0) for k in ["Indian", "Any_Other_Asian_Background", "White_And_Asian", "Unknown"]); return int(tot - nonw), int(tot)
rows = []
def add(c, role, plat, d, sexinf=None):
    fem = d.sex.str.lower() == "female"; pa = stats.ttest_ind(d.age[d.AD == 1], d.age[d.AD == 0]).pvalue; ps = stats.chi2_contingency(pd.crosstab(d.AD, fem), correction=True)[1]
    eth_s = "not reported" if c == "GSE140829" else "White/European {}/{}".format(*white(c))
    for g, lab in [(1, "AD"), (0, "CN")]:
        x = d[d.AD == g]
        r = {"Cohort": c, "Role": role, "Platform": plat, "Group": lab, "n": len(x), "Age, mean ± SD": f"{x.age.mean():.1f} ± {x.age.std():.1f}",
             "Female (annotated), n (%)": f"{int((x.sex.str.lower() == 'female').sum())} ({100 * (x.sex.str.lower() == 'female').mean():.1f}%)",
             "Female (expression-inferred), n (%)": "" if sexinf is None else f"{int((x[sexinf] == 'female').sum())} ({100 * (x[sexinf] == 'female').mean():.1f}%)",
             "APOE ε4 carriers": f"{int(x.apoe.astype(str).str.contains('4').sum())} ({100 * x.apoe.astype(str).str.contains('4').mean():.1f}%)" if c == "GSE140829" else "not available",
             "Ethnicity (whole cohort)": eth_s if g == 1 else "", "p (age; sex)": "" if g == 1 else f"{pa:.2g}; {ps:.2f}"}
        rows.append(r)
d = pd.read_csv("data/GSE140829_z.csv.gz", usecols=lambda x: x in ["id", "AD", "age", "sex", "apoe"]); add("GSE140829", "Discovery", "HumanHT-12 v4", d)
for c, plat in [("GSE63061", "HumanHT-12 v4"), ("GSE63060", "HumanHT-12 v3")]:
    w = pd.read_csv(f"data/{c}_inclSexDiscordant_z.csv.gz", usecols=lambda x: x in ["id", "AD", "age", "sex", "sex_inferred"]); add(c, "Validation", plat, w, "sex_inferred")
T1 = pd.DataFrame(rows); T1.to_csv(R_ + "Table1_cohorts.csv", index=False)
GT = pd.read_csv(RV + "AH_specificity_group_age_sex.csv")
GT.assign(age=lambda x: x.age_mean.round(1).astype(str) + " ± " + x.age_sd.round(1).astype(str), female=lambda x: x.female_pct.round(1).astype(str) + "%")[["group", "n", "age", "female"]] \
  .rename(columns={"group": "Specificity group", "age": "Age, mean ± SD", "female": "Female"}).to_csv(R_ + "Table1b_specificity_groups.csv", index=False)
# signature genes
co = pd.read_csv(R_ + "frozen_model_coefficients.csv"); G = co.term[1:].tolist(); beta = co.set_index("term").coefficient
nm = json.load(open(R_ + "signature_gene_names_NCBI.json")); TT = pd.read_csv(RV + "AK_tstats_for_enrichment.csv").set_index("gene")
fr = pd.read_csv(R_ + "discovery_selection_frequency.csv").set_index("gene").selection_frequency
gc = pd.read_csv(RV + "N_signature_gene_celltype_correlations.csv", index_col=0)
lab = {"Neutrophil": "neutrophil", "Monocyte": "monocyte", "Bcell": "B cell", "Tcell": "T cell", "NK": "NK cell", "Platelet_erythroid": "platelet/erythroid"}
pid = pd.read_csv(R_ + "signature_probe_identity.csv").set_index("gene")
def corr_lab(g):
    k = gc.loc[g].idxmax(); r = gc.loc[g].max(); return f"{lab[k]} (r = {r:.2f}{', weak' if abs(r) < 0.3 else ''})"
T = pd.DataFrame({"Gene": G, "Approved name (NCBI Gene / HGNC)": [nm[g][1] for g in G],
    "t discovery (age, sex, batch)": [round(TT.loc[g, "disc_batch_adj"], 2) for g in G], "t GSE63061 (age, sex)": [round(TT.loc[g, "GSE63061"], 2) for g in G], "t GSE63060 (age, sex)": [round(TT.loc[g, "GSE63060"], 2) for g in G],
    "CV selection %": [int(round(100 * fr.get(g, 0))) for g in G], "Strongest positive cell-type correlate": [corr_lab(g) for g in G],
    "Probe GSE140829 (nuID)": [pid.loc[g, "discovery_probe_nuID"] for g in G], "Probe GSE63061": [pid.loc[g, "GSE63061_primary_probe"] for g in G],
    "Probe GSE63060": [pid.loc[g, "GSE63060_primary_probe"] for g in G],
    "Same sequence as discovery probe (61/60)": [("yes" if pid.loc[g, "GSE63061_primary_is_same_sequence"] else "no") + "/" + ("yes" if pid.loc[g, "GSE63060_primary_is_same_sequence"] else "no") for g in G]})
T.to_csv(R_ + "Table2_signature_genes.csv", index=False)
# agreement (primary analysis sets): Spearman correlation of t-statistics; direction concordance of signature genes vs a null of
# 5,000 random 17-gene sets drawn from the non-signature genes whose discovery |t| is at least that of the weakest signature gene
rng = np.random.default_rng(20261029); tmin = TT.loc[G, "disc_age_sex"].abs().min(); pool = [g for g in TT.index[TT.disc_age_sex.abs() >= tmin] if g not in G]
ag = []
for c in ["GSE63061", "GSE63060"]:
    rho, prho = stats.spearmanr(TT.disc_age_sex, TT[c]); conc = int((np.sign(TT.loc[G, "disc_age_sex"]) == np.sign(TT.loc[G, c])).sum())
    null = np.array([(np.sign(TT.loc[s, "disc_age_sex"]) == np.sign(TT.loc[s, c])).sum() for s in [rng.choice(pool, 17, replace=False) for _ in range(5000)]])
    ag.append(dict(cohort=c, spearman_rho=rho, rho_p=prho, concordant=conc, sign_test_p=stats.binomtest(conc, 17, 0.5, alternative="greater").pvalue,
                   pool_size=len(pool), abs_t_min=tmin, null_mean=null.mean(), empirical_p=((null >= conc).sum() + 1) / 5001))
pd.DataFrame(ag).to_csv(RV + "AL_agreement_primary.csv", index=False); print(pd.DataFrame(ag))
ck = [{"file": "raw/" + os.path.basename(f), "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest(), "bytes": os.path.getsize(f)} for f in sorted(glob.glob("raw/*"))]
for f in [R_ + "frozen_model_coefficients.csv", "data/genesets/Li2014_BTM.gmt"]:
    if os.path.exists(f): ck.append({"file": f, "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest(), "bytes": os.path.getsize(f)})
pd.DataFrame(ck).to_csv(R_ + "checksums_sha256.csv", index=False)
pd.Series({"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "statsmodels": statsmodels.__version__, "matplotlib": matplotlib.__version__}, name="version").to_csv(R_ + "python_versions.csv")
print(T1.to_string())
