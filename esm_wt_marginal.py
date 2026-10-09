"""DOC-2-073 scoring code (frozen at lock-1). Wild-type-marginal ESM-2 scoring, fp32, CPU.
Per AMENDMENT-2/3. Usage:
  python esm_wt_marginal.py --measure-ram          # data-independent RAM check, dummy sequence
"""
import argparse, json, os, resource, sys, time
import numpy as np, torch
from transformers import AutoTokenizer, EsmForMaskedLM

AA = "ACDEFGHIKLMNPQRSTVWY"
WIN, STRIDE = 1022, 511
MODELS = {"650M": "facebook/esm2_t33_650M_UR50D", "150M": "facebook/esm2_t30_150M_UR50D"}

def load(size):
    tok = AutoTokenizer.from_pretrained(MODELS[size])
    model = EsmForMaskedLM.from_pretrained(MODELS[size], torch_dtype=torch.float32).eval()
    return tok, model

@torch.no_grad()
def wt_logprobs(seq, tok, model):
    """Return (L x 20) log p(aa | unmasked wt window), columns ordered as AA.
    >1022 aa: windows of 1022, stride 511; each position uses the window where it is
    farthest from a window edge (ties: earlier window)."""
    L = len(seq)
    starts = [0] if L <= WIN else list(range(0, L - WIN + 1, STRIDE))
    if L > WIN and starts[-1] + WIN < L:
        starts.append(L - WIN)
    ids = [tok.convert_tokens_to_ids(a) for a in AA]
    out = np.zeros((L, 20), dtype=np.float64)
    best = np.full(L, -1.0)
    for s in starts:
        sub = seq[s:s + WIN]
        enc = tok(sub, return_tensors="pt")
        lp = torch.log_softmax(model(**enc).logits[0, 1:-1].float(), dim=-1)[:, ids].numpy()
        n = len(sub)
        for i in range(n):
            d = min(i, n - 1 - i)
            if d > best[s + i]:
                best[s + i] = d
                out[s + i] = lp[i]
    return out

def wtll(lp, seq):
    """Mean over ALL residues of the full sequence of log p(wt aa_i), fixed per protein."""
    return float(np.mean([lp[i, AA.index(a)] for i, a in enumerate(seq)]))

def score(lp, pos1, wt, mut):
    """wt-marginal: log p(mut) - log p(wt) at 1-based position. Single substitutions only."""
    i = pos1 - 1
    return float(lp[i, AA.index(mut)] - lp[i, AA.index(wt)])

def measure(size):
    tok, model = load(size)
    rng = np.random.default_rng(0)
    seq = "".join(rng.choice(list(AA), WIN))
    wt_logprobs(seq, tok, model)
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    mem = [l for l in open("/proc/meminfo") if l.startswith("MemTotal")][0].split()[1]
    return {"model": size, "peak_rss_gb": round(peak_gb, 2), "mem_total_gb": round(int(mem) / 1024 / 1024, 2),
            "cuda": torch.cuda.is_available(), "torch": torch.__version__}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--measure-ram", action="store_true")
    a = ap.parse_args()
    if a.measure_ram:
        print("MEASURE_JSON " + json.dumps(measure("650M")))
