"""8_duplicates_ethnicity.py — identify GSE140830 arrays that duplicate any (pre-QC) GSE140829 array, using fixed-reference
scaled data jointly centred across both series; tabulate AddNeuroMed ethnicity annotations. Reuses helpers from 0_download_qc.py."""
import os; os.makedirs('results/sensitivity', exist_ok=True)
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

# ---------------- 8: duplicate check vs ALL pre-QC GSE140829 arrays; AddNeuroMed ethnicity ----------------
universe=pd.read_csv("results/discovery_gene_means.csv").gene.tolist()
ref=pd.read_csv("results/fixed_reference_parameters.csv")
ref_q=ref.ref_quantile.values; mu=ref.set_index("gene").ref_mean; sd=ref.set_index("gene").ref_sd
def refz(G):
    M=G.reindex(universe).values; R=np.empty_like(M)
    for j in range(M.shape[1]):
        r=pd.Series(M[:,j]).rank(method="average").values-1; R[:,j]=np.interp(r,np.arange(len(ref_q)),ref_q)
    return pd.DataFrame(R,index=universe,columns=G.columns).sub(mu,axis=0).div(sd,axis=0)
sym=annot_soft(fetch("GPL15988_soft"))
meta,X=read_matrix(fetch("GSE140829")); ch=chars(meta); A=refz(to_gene(X,sym)); dxA=ch["diagnosis"]
meta,X=read_matrix(fetch("GSE140830")); ch3=chars(meta); Bm=refz(to_gene(X,sym)); dxB=ch3["diagnosis"]
J=pd.concat([A,Bm],axis=1); J=J.sub(J.mean(axis=1),axis=0)          # joint centring across both series
a=J.iloc[:,:A.shape[1]].values; b=J.iloc[:,A.shape[1]:].values
a=a/np.linalg.norm(a,axis=0); b=b/np.linalg.norm(b,axis=0); C=b.T@a
best=C.max(1); arg=C.argmax(1)
out=pd.DataFrame({"id":Bm.columns,"dx":dxB.reindex(Bm.columns).values,"best_r":best.round(4),"match_GSE140829":A.columns[arg],"match_dx":dxA.reindex(A.columns[arg]).values})
out.to_csv("results/sensitivity/P_duplicates_vs_all_preQC_GSE140829.csv",index=False)
q=pd.read_csv("data/GSE140830_fixedref.csv.gz",usecols=["id","dx"]).merge(out[["id","best_r","match_GSE140829","match_dx"]],on="id",how="left")
q["dup_of"]=np.where(q.best_r>0.9,q.match_GSE140829,""); q["dup_r"]=q.best_r; q["dup_dx"]=np.where(q.best_r>0.9,q.match_dx,"")
q[["id","dx","dup_of","dup_r","dup_dx"]].to_csv("results/GSE140830_duplicates_of_GSE140829.csv",index=False)
print(pd.crosstab(out.dx,out.best_r>0.9)); print("match dx among r>0.9:",out[out.best_r>0.9].match_dx.value_counts().to_dict())
print("best r distribution (sorted, around gap):",np.sort(best)[::-1][[0,100,200,230,238,239,240,245,260]].round(3))
eth=[]
for acc in ["GSE63061","GSE63060"]:
    meta,_=read_matrix(fetch(acc)); chA=chars(meta); ids=pd.read_csv(f"data/{acc}_inclSexDiscordant_z.csv.gz",usecols=["id"]).id   # primary analysis set (sex-discordant samples retained)
    e=chA["ethnicity"].reindex(ids).value_counts().to_dict(); eth.append(dict(cohort=acc,**e))
pd.DataFrame(eth).to_csv("results/sensitivity/P_ethnicity_AddNeuroMed.csv",index=False); print(eth)
