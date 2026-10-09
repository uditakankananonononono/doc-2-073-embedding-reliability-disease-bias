"""DOC-2-073 step 1: ESM-2 650M wt-marginal scoring (fp32, CPU) of single substitutions. Frozen at lock-1 + AMENDMENT-6.
Usage: python run_scoring.py  -> per_assay.csv (one row per assay). No analysis happens here."""
import io, os, sys, zipfile, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from esm_wt_marginal import load, wt_logprobs, wtll, score, AA

ref = pd.read_csv("DMS_substitutions.csv")
excl = set(pd.read_csv("exclusions.tsv", sep="\t").DMS_id)
expo = pd.read_csv("exposures.tsv", sep="\t")
expo = expo[expo.annotation_score.notna()]
keep = [i for i in ref.DMS_id if i not in excl and i in set(expo.DMS_id)]
print("assays kept", len(keep), flush=True)
tok, model = load("650M")
cache = {}
zf = zipfile.ZipFile("DMS_ProteinGym_substitutions.zip")
names = {os.path.basename(n): n for n in zf.namelist() if n.endswith(".csv")}
rows = []
t0 = time.time()
for k, i in enumerate(keep):
    r = ref[ref.DMS_id == i].iloc[0]
    seq = r.target_seq
    if seq not in cache:
        lp = wt_logprobs(seq, tok, model)
        cache[seq] = (lp, wtll(lp, seq))
    lp, wl = cache[seq]
    d = pd.read_csv(io.BytesIO(zf.read(names[r.DMS_filename])))
    d = d[~d.mutant.str.contains(":")]
    s, y, bad = [], [], 0
    for m, v in zip(d.mutant, d.DMS_score):
        wt, pos, mu = m[0], int(m[1:-1]), m[-1]
        if pos < 1 or pos > len(seq) or seq[pos - 1] != wt or mu not in AA or wt not in AA:
            bad += 1; continue
        s.append(score(lp, pos, wt, mu)); y.append(v)
    ok = len(s) >= 50
    rho = spearmanr(s, y)[0] if ok else np.nan
    rows.append(dict(DMS_id=i, n_scored=len(s), n_mismatch=bad, rho=rho, wtll=wl, seq_len=len(seq)))
    print(k, i, len(s), round(time.time() - t0), flush=True)
    pd.DataFrame(rows).to_csv("per_assay.csv", index=False)
print("SCORING_DONE", len(rows), flush=True)
