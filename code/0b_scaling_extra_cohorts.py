"""0b_scaling_extra_cohorts.py — alternative scalings (control-referenced; fixed-reference single-sample) for discovery and
validation cohorts, plus GSE140829 MCI and GSE140830 (FTD spectrum) for the specificity analysis. Reuses helper functions from 0_download_qc.py."""
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

# ---------------- 0b: alternative scalings + additional groups (MCI, FTD spectrum) ----------------
U["GSE140830"]="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140830/matrix/GSE140830_series_matrix.txt.gz"
universe=pd.read_csv("results/discovery_gene_means.csv").gene.tolist()
def qc_simple(name,G,ph):
    ysc=G.loc[[g for g in ["RPS4Y1","DDX3Y","KDM5D","EIF1AY"] if g in G.index]].mean(); xs=G.loc["XIST"]
    d=(ysc-ysc.mean())/ysc.std()-(xs-xs.mean())/xs.std(); inf=np.where(d>0,"male","female")
    disc=inf!=ph["sex"].str.lower().values
    C=np.corrcoef(G.values.T); np.fill_diagonal(C,np.nan); mc=np.nanmean(C,axis=1)
    med=np.median(mc); mad=np.median(np.abs(mc-med))*1.4826; out=mc<med-4*mad
    keep=~(disc|out); log.append(f"{name}: sex-discordant {disc.sum()}, outliers {out.sum()}, kept {keep.sum()}/{len(keep)}")
    return G.loc[:,keep], ph[keep]
sym15988=annot_soft(fetch("GPL15988_soft"))
meta,X=read_matrix(fetch("GSE140829")); ch=chars(meta)
ph=pd.DataFrame({"id":ch.index,"dx":ch["diagnosis"],"age":pd.to_numeric(ch["age_at_draw"],errors="coerce"),"sex":ch["sex"],"batch":ch["batch"]},index=ch.index)
Gall=to_gene(X[ph.index],sym15988)
keep_adcn=pd.read_csv("data/GSE140829_z.csv.gz",usecols=["id"]).id.tolist()
mci_G,mci_ph=qc_simple("GSE140829 MCI",Gall.loc[:,ph.index[ph.dx=="MCI"]],ph[ph.dx=="MCI"])
cohorts={"GSE140829":(Gall[keep_adcn],ph.loc[keep_adcn])}
cohorts["GSE140829_MCI"]=(mci_G,mci_ph)
for a,gpl in [("GSE63061","GPL10558"),("GSE63060","GPL6947")]:
    meta,X=read_matrix(fetch(a)); ch=chars(meta)
    p2=pd.DataFrame({"id":ch.index,"dx":ch["status"],"age":pd.to_numeric(ch["age"],errors="coerce"),"sex":ch["gender"]},index=ch.index)
    ids=pd.read_csv(f"data/{a}_z.csv.gz",usecols=["id"]).id.tolist()
    Gv=to_gene(X,annot_annot(fetch(gpl))); cohorts[a]=(Gv[ids],p2.loc[ids])
    # analysis set retaining sex-discordant samples (outliers/duplicates removed); used as the primary external set
    ids_r=pd.read_csv(f"data/{a}_inclSexDiscordant_z.csv.gz",usecols=["id"]).id.tolist(); cohorts[a+"_retained"]=(Gv[ids_r],p2.loc[ids_r])
    Gv[ids_r].T.reset_index().rename(columns={"index":"id"}).to_csv(f"data/{a}_retained_log2_allgenes.csv.gz",index=False,compression="gzip")
meta,X=read_matrix(fetch("GSE140830")); ch=chars(meta)
log.append("GSE140830 characteristics: "+"; ".join(f"{c}: {dict(ch[c].value_counts().head(12))}" for c in ch.columns if ch[c].nunique()<15))
dxcol=[c for c in ch.columns if "diag" in c][0]; agecol=[c for c in ch.columns if "age" in c][0]
p3=pd.DataFrame({"id":ch.index,"dx":ch[dxcol],"age":pd.to_numeric(ch[agecol],errors="coerce"),"sex":ch["sex"],"batch":ch.get("batch")},index=ch.index)
cohorts["GSE140830"]=qc_simple("GSE140830",to_gene(X,sym15988),p3)
# reference distribution from discovery AD/CN
D=cohorts["GSE140829"][0].loc[universe]
ref_q=np.sort(D.values,axis=0).mean(axis=1)
def refnorm(G):
    M=G.loc[universe].values; R=np.empty_like(M)
    for j in range(M.shape[1]):
        r=pd.Series(M[:,j]).rank(method="average").values-1
        R[:,j]=np.interp(r,np.arange(len(ref_q)),ref_q)
    return pd.DataFrame(R,index=universe,columns=G.columns)
Dq=refnorm(cohorts["GSE140829"][0]); mu=Dq.mean(axis=1); sd=Dq.std(axis=1)
pd.DataFrame({"gene":universe,"ref_quantile":ref_q,"ref_mean":mu.values,"ref_sd":sd.values}).to_csv("results/fixed_reference_parameters.csv",index=False)
def save(name,ph,Z,tag):
    out=pd.concat([ph[[c for c in ["id","dx","age","sex","batch"] if c in ph]].reset_index(drop=True),Z.T.reset_index(drop=True)],axis=1)
    out.to_csv(f"data/{name}_{tag}.csv.gz",index=False,compression="gzip")
for name,(G,ph) in cohorts.items():
    G=G.reindex(universe)
    Zf=refnorm(G.fillna(G.mean(axis=1).mean())).sub(mu,axis=0).div(sd,axis=0); save(name,ph,Zf,"fixedref")
    ctl=ph.dx.astype(str).str.lower().isin(["control","ctl"]).values
    if ctl.sum()>10:
        Zc=G.sub(G.loc[:,ctl].mean(axis=1),axis=0).div(G.loc[:,ctl].std(axis=1),axis=0); save(name,ph,Zc,"ctlref")
        Zz=G.sub(G.mean(axis=1),axis=0).div(G.std(axis=1),axis=0); save(name,ph,Zz,"cohortz")
    log.append(f"{name}: {G.shape[1]} samples, missing universe genes {int(G.isna().all(axis=1).sum())}")
# duplicates between GSE140830 and GSE140829 are identified in 8_duplicates_ethnicity.py (joint centring, all pre-QC arrays)
open("data/qc_summary_0b.txt","w").write("\n".join(log)); print("\n".join(log))
