"""DOC-2-073 data acquisition. Reads ONLY metadata (reference file + UniProt). Never opens DMS score files.
Writes DATA_HASHES.tsv, exclusions.tsv, exposures.tsv."""
import hashlib, io, sys, time, urllib.request, urllib.parse, csv, zipfile
import pandas as pd
Z = "https://zenodo.org/records/15293562/files/%s?download=1"
EXPECT_MD5 = {"DMS_substitutions.csv": "c434631737013fceb56efc98056151e0",
              "DMS_ProteinGym_substitutions.zip": "ca1a4d46941ef33cc972245347118c7d"}
def get(url, dest):
    urllib.request.urlretrieve(url, dest)
def h(path, algo):
    m = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): m.update(b)
    return m.hexdigest()
rows = []
for fn in EXPECT_MD5:
    import os
    if not os.path.exists(fn): get(Z % fn, fn)
    md5, sha = h(fn, "md5"), h(fn, "sha256")
    rows.append((fn, "zenodo.org/records/15293562", md5, md5 == EXPECT_MD5[fn], sha))
ref = pd.read_csv("DMS_substitutions.csv")
ids = sorted(ref.UniProt_ID.unique())
recs, rel = {}, "unknown"
for i in range(0, len(ids), 40):
    q = "(" + " OR ".join("id:" + x for x in ids[i:i + 40]) + ")"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "format": "tsv", "size": 500, "fields": "id,accession,annotation_score,lit_pubmed_id,sequence,reviewed"})
    with urllib.request.urlopen(url) as r:
        rel = r.headers.get("X-UniProt-Release", rel)
        txt = r.read().decode()
    for line in txt.strip().split("\n")[1:]:
        p = line.split("\t")
        if len(p) < 6:
            print("SHORTROW", len(p), line[:80]); p = p + [""] * (6 - len(p))
        recs[p[0]] = p
raw = "\n".join("\t".join(recs[k]) for k in sorted(recs))
open("uniprot_raw.tsv", "w").write(raw)
rows.append(("uniprot_raw.tsv", "rest.uniprot.org release " + rel, "", "", hashlib.sha256(raw.encode()).hexdigest()))
with open("DATA_HASHES.tsv", "w") as f:
    f.write("file\tsource\tmd5\tmd5_matches_zenodo\tsha256\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
ex, expo = [], []
for _, r in ref.iterrows():
    why = []
    if r.DMS_number_single_mutants < 50: why.append("lt50_single_substitutions")
    if r.UniProt_ID not in recs: why.append("uniprot_entry_not_found")
    if why: ex.append((r.DMS_id, ";".join(why)))
    p = recs.get(r.UniProt_ID)
    if p:
        n = len([x for x in p[3].split(";") if x.strip()]) if p[3] else 0
        expo.append((r.DMS_id, r.UniProt_ID, p[1], p[2], n, p[5], r.year, r.taxon, r.seq_len, r.includes_multiple_mutants))
open("exclusions.tsv", "w").write("DMS_id\treason\n" + "".join("%s\t%s\n" % e for e in ex))
open("exposures.tsv", "w").write("DMS_id\tUniProt_ID\taccession\tannotation_score\tn_pubmed\treviewed\tyear\ttaxon\tseq_len\thas_multi\n" +
    "".join("\t".join(map(str, e)) + "\n" for e in expo))
print("ACQ_DONE release=%s assays=%d proteins=%d found=%d excluded=%d" % (rel, len(ref), len(ids), len(recs), len(ex)))
for r in rows: print("HASH", *r)
