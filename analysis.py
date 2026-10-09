"""DOC-2-073 step 2: frozen analysis G1-G3 (AMENDMENT-6). Run once after per_assay.csv is complete."""
import numpy as np, pandas as pd
from scipy.stats import rankdata
rng = np.random.default_rng(12345)
a = pd.read_csv("per_assay.csv"); e = pd.read_csv("exposures.tsv", sep="\t").drop(columns=["seq_len"])
d = a.merge(e, on="DMS_id"); d = d[d.rho.notna() & (d.n_scored >= 50)].copy()
d["x"] = np.log1p(d.n_pubmed); d["loglen"] = np.log(d.seq_len)
disc, ext = d[d.year <= 2021].copy(), d[d.year >= 2022].copy()
TAX = ["Human", "Eukaryote", "Prokaryote", "Virus"]

def resid(v, X):
    X1 = np.column_stack([np.ones(len(X)), X]); b = np.linalg.lstsq(X1, v, rcond=None)[0]; return v - X1 @ b
def partial(df):
    X = np.column_stack([rankdata(df.loglen)] + [(df.taxon == t).astype(float) for t in TAX[1:]])
    return np.corrcoef(resid(rankdata(df.x), X), resid(rankdata(df.rho), X))[0, 1]
def boot(df, n=10000):
    prots = df.UniProt_ID.unique(); grp = {p: df[df.UniProt_ID == p] for p in prots}; out = []
    for _ in range(n):
        s = pd.concat([grp[p] for p in rng.choice(prots, len(prots))]); 
        if s.x.nunique() > 1: out.append(partial(s))
    return np.percentile(out, [2.5, 97.5])

res = {}
g1 = partial(disc); ci1 = boot(disc); res["G1"] = dict(n=len(disc), partial_rho=g1, ci=list(ci1), pass_=bool(g1 >= 0.20 and ci1[0] > 0))
if len(ext) >= 40:
    g2 = partial(ext); ci2 = boot(ext)
    res["G2"] = dict(n=len(ext), partial_rho=g2, ci=list(ci2), pass_=bool(np.sign(g2) == np.sign(g1) and g2 >= 0.15 and ci2[0] > 0))
else:
    res["G2"] = dict(n=len(ext), label="INCONCLUSIVE-UNDERPOWERED", pass_=False)

# G3: leave-one-taxon-group-out on discovery, nonviral folds in the gate; virus fold reported only.
def feats(df, with_expo):
    cols = [df.loglen.values, df.wtll.values]
    if with_expo: cols += [df.x.values, df.annotation_score.values]
    return np.column_stack(cols)
def ridge(Xtr, ytr, Xte, lam=1.0):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9; A = (Xtr - mu) / sd; B = (Xte - mu) / sd
    w = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ (ytr - ytr.mean())); return B @ w + ytr.mean()
def lto_r2(df, with_expo, folds):
    sse = sst = 0.0
    for t in folds:
        tr, te = df[df.taxon != t], df[df.taxon == t]
        if len(te) < 5 or len(tr) < 10: continue
        p = ridge(feats(tr, with_expo), tr.rho.values, feats(te, with_expo))
        sse += ((te.rho.values - p) ** 2).sum(); sst += ((te.rho.values - tr.rho.mean()) ** 2).sum()
    return 1 - sse / sst
folds = ["Human", "Eukaryote", "Prokaryote"]
base, full = lto_r2(disc, False, folds), lto_r2(disc, True, folds); delta = full - base
prot = disc.drop_duplicates("UniProt_ID")[["UniProt_ID", "x", "annotation_score"]].reset_index(drop=True); perm = []
for _ in range(2000):
    sh = prot.sample(frac=1, random_state=int(rng.integers(1e9))).reset_index(drop=True)
    m = dict(zip(prot.UniProt_ID, zip(sh.x, sh.annotation_score))); dd = disc.copy()
    dd["x"] = [m[u][0] for u in dd.UniProt_ID]; dd["annotation_score"] = [m[u][1] for u in dd.UniProt_ID]
    perm.append(lto_r2(dd, True, folds) - base)
pv = (1 + sum(p >= delta for p in perm)) / (1 + len(perm))
res["G3"] = dict(base_r2=base, full_r2=full, delta=delta, perm_p=pv, pass_=bool(delta >= 0.03 and pv < 0.05),
                 virus_fold_exploratory=dict(base=lto_r2(disc, False, ["Virus"]), full=lto_r2(disc, True, ["Virus"])))
res["LABEL"] = "POSITIVE" if all(res[g]["pass_"] for g in ["G1", "G2", "G3"]) else ("INCONCLUSIVE-UNDERPOWERED" if res["G2"].get("label") else "HONEST NEGATIVE")
import json; print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))
