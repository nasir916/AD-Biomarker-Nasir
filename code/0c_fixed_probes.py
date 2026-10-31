
"""0c_fixed_probes.py — probe-level transfer checks for the signature genes.
(1) decode the lumi nuIDs of GSE140829 (GPL15988) to probe sequences and map each to the ILMN probe with the
identical 50-nt sequence on GPL10558 (GSE63061) and GPL6947 (GSE63060);
(2) build validation matrices in which every signature gene is measured by the probe with the SAME sequence as the
discovery probe (fixed-probe sensitivity), and (3) by the average of all probes annotated to the gene (probe-averaging).
Matrices are z-scored within cohort over the QC-passed samples, as in the primary analysis."""
import gzip, collections, numpy as np, pandas as pd
CODE=[chr(c) for c in range(65,91)]+[chr(c) for c in range(97,123)]+list("0123456789")+["_","."]
IX={c:i for i,c in enumerate(CODE)}
def nuid2seq(n):                      # lumi nuID: 1 check/padding character + 6-bit codes of nucleotide triplets
    return "".join("ACGT"[IX[c]//16]+"ACGT"[(IX[c]//4)%4]+"ACGT"[IX[c]%4] for c in n[1:])
def readplat(p):
    rows=[];hdr=None;on=False
    with gzip.open(p,"rt",errors="ignore") as f:
        for l in f:
            if l.startswith("!platform_table_begin"): on=True; continue
            if l.startswith("!platform_table_end"): break
            if on:
                x=l.rstrip("\n").split("\t")
                if hdr is None: hdr=x
                else: rows.append(x)
    return pd.DataFrame(rows,columns=hdr)
def read_matrix(path):
    rows=[];hdr=None;on=False
    with gzip.open(path,"rt",errors="ignore") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"): on=True; continue
            if l.startswith("!series_matrix_table_end"): break
            if on:
                x=[v.strip('"') for v in l.rstrip("\n").split("\t")]
                if hdr is None: hdr=x
                else: rows.append(x)
    X=pd.DataFrame(rows,columns=hdr).set_index(hdr[0]).apply(pd.to_numeric,errors="coerce")
    return X                                   # values used as deposited, as in the primary analysis
P={"GSE63061":readplat("raw/GPL10558_platform_full.txt.gz"),"GSE63060":readplat("raw/GPL6947_platform_full.txt.gz")}
co=pd.read_csv("results/frozen_model_coefficients.csv"); genes=[g for g in co.iloc[:,0] if g!="Intercept"]
pm=pd.read_csv("results/probe_map_by_cohort.csv").set_index("gene")
# genome-wide sequence map
Xd_ids=read_matrix("raw/GSE140829.gz").index
seqmap=pd.DataFrame({"nuID":Xd_ids,"sequence":[nuid2seq(n)[:50] for n in Xd_ids]})
for a,pl in P.items():
    d={s.upper():i for i,s in zip(pl.ID,pl.SEQUENCE) if s}; seqmap["ILMN_"+a]=seqmap.sequence.map(d)
seqmap.to_csv("results/probe_sequence_map_GSE140829.csv.gz",index=False,compression="gzip")
summ=[f"GSE140829 probes decoded: {len(seqmap)}; identical-sequence probe on GPL10558: {seqmap.ILMN_GSE63061.notna().sum()}; on GPL6947: {seqmap.ILMN_GSE63060.notna().sum()}"]
rows=[]
for g in genes:
    nu=pm.loc[g,"probe_GSE140829"]; r=seqmap.set_index("nuID").loc[nu]
    rows.append({"gene":g,"discovery_probe_nuID":nu,"probe_sequence":r.sequence,
      "GSE63061_primary_probe":pm.loc[g,"probe_GSE63061"],"GSE63061_same_sequence_probe":r.ILMN_GSE63061,
      "GSE63060_primary_probe":pm.loc[g,"probe_GSE63060"],"GSE63060_same_sequence_probe":r.ILMN_GSE63060})
T=pd.DataFrame(rows)
for a in P: T[a+"_primary_is_same_sequence"]=T[a+"_primary_probe"]==T[a+"_same_sequence_probe"]
T.to_csv("results/signature_probe_identity.csv",index=False)
for a,pl in P.items():
    X=read_matrix(f"raw/{a}.gz"); q=pd.read_csv(f"data/{a}_z.csv.gz",usecols=["id"]).id
    X=X[q]
    import io
    gp={"GSE63061":"GPL10558","GSE63060":"GPL6947"}[a]
    with gzip.open(f"raw/{gp}.gz","rt",errors="ignore") as f: tx=f.read()
    an=pd.read_csv(io.StringIO(tx.split("!platform_table_begin\n")[1].split("!platform_table_end")[0]),sep="\t",low_memory=False)
    sym=an.set_index("ID")["Gene symbol"].reindex(X.index)   # same annotation as the primary probe-to-gene mapping
    fx={}; av={}; npr={}
    for _,r in T.iterrows():
        g=r.gene; pr=r[a+"_same_sequence_probe"] if isinstance(r[a+"_same_sequence_probe"],str) else r[a+"_primary_probe"]
        fx[g]=X.loc[pr]
        allp=sym.index[sym==g]; allp=[p for p in allp if p in X.index and X.loc[p].notna().all()]
        npr[g]=len(allp); av[g]=X.loc[allp].mean()
    for tag,D in [("fixedprobe",fx),("probemean",av)]:
        Z=pd.DataFrame(D); Z=(Z-Z.mean())/Z.std()
        Z.insert(0,"id",Z.index); Z.to_csv(f"data/{a}_{tag}_z.csv.gz",index=False,compression="gzip")
    T[a+"_n_probes_for_gene"]=T.gene.map(npr)
    summ.append(f"{a}: primary probe identical in sequence to discovery probe for {int(T[a+'_primary_is_same_sequence'].sum())}/{len(T)} genes; same-sequence probe available for {int(T[a+'_same_sequence_probe'].notna().sum())}/{len(T)}")
T.to_csv("results/signature_probe_identity.csv",index=False)
open("results/probe_identity_summary.txt","w").write("\n".join(summ)); print("\n".join(summ))
