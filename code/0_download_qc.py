
"""0_download_qc.py — download GSE140829 (discovery), GSE63061 and GSE63060 (validation) from GEO,
map probes to genes (highest-mean probe per gene within cohort), QC (expression-inferred sex,
array outliers, duplicates), restrict to the shared gene universe, z-score within cohort.
Cohorts: Illumina HumanHT-12 whole-blood series with >= 40 AD and >= 40 CN (see manuscript Methods 2.2)."""
import os, io, gzip, requests, collections, numpy as np, pandas as pd
RAW="raw"; os.makedirs(RAW,exist_ok=True); os.makedirs("data",exist_ok=True)
U={"GSE140829":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140829/matrix/GSE140829_series_matrix.txt.gz",
 "GSE63060":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63060/matrix/GSE63060_series_matrix.txt.gz",
 "GSE63061":"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63061/matrix/GSE63061_series_matrix.txt.gz",
 "GPL15988_soft":"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15988/soft/GPL15988_family.soft.gz",
 "GPL6947":"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz",
 "GPL10558":"https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"}
def fetch(k):
    p=f"{RAW}/{k}.gz"
    if not os.path.exists(p):
        r=requests.get(U[k],timeout=1800); r.raise_for_status(); open(p,"wb").write(r.content)
    return p
def read_matrix(path):
    meta={}
    with gzip.open(path,"rt",errors="ignore") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): break
            if l.startswith("!Sample_"):
                k=l.split("\t")[0][8:]; v=[x.strip().strip('"') for x in l.rstrip("\n").split("\t")[1:]]
                meta.setdefault(k,[]).append(v)
        X=pd.read_csv(f,sep="\t",index_col=0,comment="!",low_memory=False)
    return meta, X.apply(pd.to_numeric,errors="coerce")
def chars(meta):
    out={}
    for row in meta["characteristics_ch1"]:
        for i,v in enumerate(row):
            if ": " in v:
                k,val=v.split(": ",1); out.setdefault(k.strip().lower(),{})[i]=val.strip()
    n=len(meta["geo_accession"][0])
    return pd.DataFrame({k:[d.get(i) for i in range(n)] for k,d in out.items()},index=meta["geo_accession"][0])
def annot_soft(p):
    with gzip.open(p,"rt",errors="ignore") as f: t=f.read()
    s=t.split("!platform_table_begin\n")[1].split("!platform_table_end")[0]
    g=pd.read_csv(io.StringIO(s),sep="\t",low_memory=False)
    return g.set_index("ID")["GENE SYMBOL"]
def annot_annot(p):
    with gzip.open(p,"rt",errors="ignore") as f: t=f.read()
    s=t.split("!platform_table_begin\n")[1].split("!platform_table_end")[0]
    g=pd.read_csv(io.StringIO(s),sep="\t",low_memory=False)
    return g.set_index("ID")["Gene symbol"]
def to_gene(X, sym):
    sym=sym.reindex(X.index)
    ok=sym.notna() & ~sym.astype(str).str.contains("///") & (sym.astype(str).str.len()>0)
    X=X[ok]; sym=sym[ok]
    mu=X.mean(axis=1)
    best=mu.groupby(sym).idxmax()
    G=X.loc[best.values]; G.index=best.index
    G.attrs["probe"]=best            # gene -> probe ID chosen (highest mean in this cohort)
    return G
log=[]
def qc(name, G, ph):
    n0=G.shape[1]
    # expression-inferred sex
    ysc=G.loc[[g for g in ["RPS4Y1","DDX3Y","KDM5D","EIF1AY"] if g in G.index]].mean()
    xs=G.loc["XIST"]; d=(ysc-ysc.mean())/ysc.std()-(xs-xs.mean())/xs.std()
    inf=np.where(d>0,"male","female")
    ann=ph["sex"].str.lower().values
    disc=inf!=ann
    log.append(f"{name}: sex-discordant {disc.sum()} (AD {int((disc&(ph.AD==1)).sum())}, CN {int((disc&(ph.AD==0)).sum())})")
    # array outliers: mean inter-array correlation
    C=np.corrcoef(G.values.T); np.fill_diagonal(C,np.nan); mc=np.nanmean(C,axis=1)
    med=np.median(mc); mad=np.median(np.abs(mc-med))*1.4826
    out=mc<med-4*mad
    log.append(f"{name}: correlation outliers {out.sum()}")
    # duplicates: correlation of gene-centred profiles
    Gc=G.sub(G.mean(axis=1),axis=0); Cc=np.corrcoef(Gc.values.T); np.fill_diagonal(Cc,0)
    iu=np.argwhere(np.triu(Cc)>0.90); dup=np.zeros(n0,bool)
    for i,j in iu: dup[j]=True
    log.append(f"{name}: duplicate pairs (centred r>0.90) {len(iu)}; raw-r>0.99 pairs {(np.triu(C,1)>0.99).sum()}")
    keep=~(disc|out|dup)
    log.append(f"{name}: kept {keep.sum()}/{n0} (AD {int(ph.AD[keep].sum())}, CN {int((ph.AD[keep]==0).sum())})")
    ph=ph.assign(sex_inferred=inf,sex_discordant=disc,outlier=out,duplicate=dup)
    return G.loc[:,keep], ph[keep], ph
# discovery
meta,X=read_matrix(fetch("GSE140829")); ch=chars(meta)
ph=pd.DataFrame({"id":ch.index,"dx":ch["diagnosis"],"age":pd.to_numeric(ch["age_at_draw"],errors="coerce"),
                 "sex":ch["sex"],"batch":ch["batch"],"apoe":ch["apoe_status"]},index=ch.index)
ph=ph[ph.dx.isin(["AD","Control"])]; ph["AD"]=(ph.dx=="AD").astype(int)
G=to_gene(X[ph.index], annot_soft(fetch("GPL15988_soft"))); PROBE_DISC=G.attrs["probe"]
disc_G, disc_ph, disc_full = qc("GSE140829", G, ph)
val={}; RAWG={}; PROBES={}
for a,gpl in [("GSE63061","GPL10558"),("GSE63060","GPL6947")]:
    meta,X=read_matrix(fetch(a)); ch=chars(meta)
    ph=pd.DataFrame({"id":ch.index,"dx":ch["status"],"age":pd.to_numeric(ch["age"],errors="coerce"),
                     "sex":ch["gender"]},index=ch.index)
    ph=ph[ph.dx.isin(["AD","CTL"])]; ph["AD"]=(ph.dx=="AD").astype(int)
    G=to_gene(X[ph.index], annot_annot(fetch(gpl)))
    val[a]=qc(a,G,ph); RAWG[a]=G; PROBES[a]=G.attrs["probe"]
universe=sorted(set(disc_G.index)&set(val["GSE63061"][0].index)&set(val["GSE63060"][0].index))
mu=disc_G.loc[universe].mean(axis=1); universe=sorted(mu[mu>=mu.quantile(0.25)].index)
log.append(f"Shared gene universe after removing lowest-quartile discovery expression: {len(universe)}")
def zsave(name,G,ph,full):
    Z=G.loc[universe]; Z=Z.sub(Z.mean(axis=1),axis=0).div(Z.std(axis=1),axis=0)
    out=pd.concat([ph[[c for c in ["id","AD","age","sex","batch","apoe"] if c in ph]].reset_index(drop=True),
                   Z.T.reset_index(drop=True)],axis=1)
    out.to_csv(f"data/{name}_z.csv.gz",index=False,compression="gzip"); full.to_csv(f"data/{name}_qc_log.csv",index=False)
    log.append(f"{name}: saved {out.shape}")
zsave("GSE140829",disc_G,disc_ph,disc_full)
os.makedirs("results",exist_ok=True)
m=disc_G.loc[universe].mean(axis=1)
pd.DataFrame({"gene":m.index,"mean_log2":m.values,"bin":pd.qcut(m,10,labels=False).values}).to_csv("results/discovery_gene_means.csv",index=False)
for a,(G,ph,full) in val.items(): zsave(a,G,ph,full)
# ---- probe actually used per gene in each cohort (reported for the signature genes) ----
pd.DataFrame({"gene":PROBE_DISC.index,"probe_GSE140829":PROBE_DISC.values}).merge(
  pd.DataFrame({"gene":PROBES["GSE63061"].index,"probe_GSE63061":PROBES["GSE63061"].values}),on="gene",how="outer").merge(
  pd.DataFrame({"gene":PROBES["GSE63060"].index,"probe_GSE63060":PROBES["GSE63060"].values}),on="gene",how="outer").to_csv("results/probe_map_by_cohort.csv",index=False)
# ---- sensitivity datasets keeping sex-discordant samples (outliers/duplicates still removed); inferred sex added ----
for a,(G,ph,full) in val.items():
    keep2=~(full.outlier.values|full.duplicate.values); Gk=RAWG[a].loc[:,keep2]; pk=full[keep2]
    Z=Gk.loc[universe]; Z=Z.sub(Z.mean(axis=1),axis=0).div(Z.std(axis=1),axis=0)
    out=pd.concat([pk[["id","AD","age","sex","sex_inferred","sex_discordant"]].reset_index(drop=True),Z.T.reset_index(drop=True)],axis=1)
    out.to_csv(f"data/{a}_inclSexDiscordant_z.csv.gz",index=False,compression="gzip")
    log.append(f"{a}: sensitivity set incl. sex-discordant {out.shape[0]} (AD {int(out.AD.sum())}, CN {int((out.AD==0).sum())})")
# ---- log2 expression of all mapped genes (analysed samples) for reference-based deconvolution ----
for name,G in [("GSE140829",disc_G),("GSE63061",val["GSE63061"][0]),("GSE63060",val["GSE63060"][0])]:
    G.T.reset_index().rename(columns={"index":"id"}).to_csv(f"data/{name}_log2_allgenes.csv.gz",index=False,compression="gzip")
# ---- duplicates across the two AddNeuroMed series (gene-centred within series, then correlated across series) ----
Za=val["GSE63061"][0].loc[universe]; Za=Za.sub(Za.mean(axis=1),axis=0); Zb=val["GSE63060"][0].loc[universe]; Zb=Zb.sub(Zb.mean(axis=1),axis=0)
A=Za.values/np.linalg.norm(Za.values,axis=0); Bm=Zb.values/np.linalg.norm(Zb.values,axis=0); Cab=Bm.T@A
pd.DataFrame({"GSE63060_id":Zb.columns,"best_match_GSE63061":Za.columns[Cab.argmax(1)],"best_r":Cab.max(1).round(4)}).to_csv("results/duplicates_GSE63060_vs_GSE63061.csv",index=False)
log.append(f"GSE63060 vs GSE63061: max cross-series gene-centred r = {Cab.max():.3f}; pairs r>0.90: {(Cab>0.9).sum()}")
open("data/qc_summary.txt","w").write("\n".join(log)); print("\n".join(log))
