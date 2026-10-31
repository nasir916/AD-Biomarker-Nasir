
"""6_site_stratified.py — site-stratified validation of the frozen model in AddNeuroMed (site = first 3 letters of GEO sample description)."""
import gzip, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
fv = pd.read_csv("results/validation_predictions.csv"); out = []
for a in ["GSE63061", "GSE63060"]:
    md = {}
    with gzip.open(f"raw/{a}.gz", "rt", errors="ignore") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): break
            if l.startswith("!Sample_geo_accession") or l.startswith("!Sample_description"):
                md[l.split("\t")[0]] = [x.strip().strip('"') for x in l.rstrip("\n").split("\t")[1:]]
    site = pd.Series([d[:3].upper() for d in md["!Sample_description"]], index=md["!Sample_geo_accession"])
    v = fv[fv.cohort == a].copy(); v["site"] = site.reindex(v.id).values
    ct = pd.crosstab(v.site, v.AD); num = den = 0; per = []
    for s, g in v.groupby("site"):
        x1 = g.lp[g.AD == 1].values; x0 = g.lp[g.AD == 0].values
        if len(x1) and len(x0): o = np.subtract.outer(x1, x0); num += (o > 0).sum() + 0.5 * (o == 0).sum(); den += o.size
        if len(x1) >= 5 and len(x0) >= 5: per.append(f"{s}: {len(x1)}/{len(x0)} AUC {(np.subtract.outer(x1, x0) > 0).mean():.2f}")
    chi = stats.chi2_contingency(ct); V = np.sqrt(chi[0] / (ct.values.sum() * (min(ct.shape) - 1)))
    v["lpz"] = (v.lp - v.lp.mean()) / v.lp.std(); vv = v[v.site.isin(ct.index[(ct > 0).all(axis=1)])]
    m = smf.logit("AD ~ lpz + C(site)", vv).fit(disp=0); ci = m.conf_int().loc["lpz"]
    out.append(dict(cohort=a, n_sites=ct.shape[0], sites_with_both=int((ct > 0).all(axis=1).sum()), cramers_V_site_AD=round(V, 3), site_AD_p=chi[1],
                    site_stratified_AUC=round(num / den, 3), n_in_stratified=len(vv), OR_per_SD_site_adj=round(np.exp(m.params.lpz), 2),
                    OR_lo=round(np.exp(ci[0]), 2), OR_hi=round(np.exp(ci[1]), 2), p=m.pvalues.lpz, per_site="; ".join(per)))
pd.DataFrame(out).to_csv("results/sensitivity/G_validation_site_stratified.csv", index=False); print(pd.DataFrame(out).T)

# ---- where do the sex-discordant GSE63060 samples come from? site by diagnosis and discordance status ----
md = {}
with gzip.open("raw/GSE63060.gz", "rt", errors="ignore") as f:
    for l in f:
        if l.startswith("!series_matrix_table_begin"): break
        if l.startswith("!Sample_geo_accession") or l.startswith("!Sample_description"):
            md[l.split("\t")[0]] = [x.strip().strip('"') for x in l.rstrip("\n").split("\t")[1:]]
site = pd.Series([d[:3].upper() for d in md["!Sample_description"]], index=md["!Sample_geo_accession"])
w = pd.read_csv("data/GSE63060_inclSexDiscordant_z.csv.gz", usecols=["id", "AD", "sex", "sex_inferred", "sex_discordant"])
w["site"] = site.reindex(w.id).values; w["group"] = np.where(w.AD == 1, "AD", np.where(w.sex_discordant, "CN sex-discordant", "CN concordant"))
tab = pd.crosstab(w.site, w.group); tab.to_csv("results/sensitivity/U_GSE63060_sex_discordant_by_site.csv")
cn = w[w.AD == 0]; chi = stats.chi2_contingency(pd.crosstab(cn.site, cn.sex_discordant))
pat = pd.crosstab(cn.sex.str.lower(), cn.sex_inferred); pat.to_csv("results/sensitivity/U_GSE63060_controls_annotated_vs_inferred_sex.csv")
print(tab); print("controls: discordance vs site chi2 p =", chi[1]); print(pat)
