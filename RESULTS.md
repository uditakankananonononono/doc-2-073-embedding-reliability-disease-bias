# DOC-2-073 results (frozen analysis, AMENDMENT-6; no post-hoc changes)

Label: HONEST NEGATIVE, CONSISTENT-WITH-NULL. No gate passed. This is "not detected, imprecise", not evidence of absence. No claim of an annotation-bias effect is made.

Independent review: scoped pass as an honest negative. The reviewer re-ran analysis.py unchanged on the committed CSVs and reproduced every number in results.json to <1e-9 (fixed seeds). This file was then amended in wording only (four items below); the label did not change.

## Results
Scoring: ESM-2 650M wild-type-marginal, fp32, Colab CPU, 204 assays. Files: per_assay.csv (rho per assay, WT log-likelihood, length), results.json.
- G1 (discovery, year <= 2021, n=98 assays after the n_scored >= 50 filter): partial rho -0.068, 95% protein-cluster bootstrap CI [-0.251, 0.140]. Fail. The CI excludes rho >= +0.20.
- G2 (external, year >= 2022, n=106): partial rho +0.101, CI [-0.137, 0.318]. Fail. The CI still contains the +0.15 pass threshold, and the unadjusted external rho is +0.216 (per the independent reviewer's reproduction). The estimate is imprecise; an effect of the size the protocol cared about is not ruled out. Sign differs from G1; both CIs span 0.
- G3 (leave-one-taxon-out on Human/Eukaryote/Prokaryote): baseline R2 -0.210, full R2 -0.252, delta -0.042, permutation p 0.544. Fail. The baseline R2 is negative, so the delta compares two poor models and says little. Virus fold (exploratory only): base 0.280, full 0.148.

## Deviations from PROTOCOL.md
1. Scoring is wild-type-marginal, not masked-marginal. This is a different estimand.
2. Compute is CPU fp32, not T4 fp16 (Colab T4 was refused for usage limits).
3. Clusters are UniProt_ID (protein-level), not InterPro-family clusters. Weaker control for relatedness (AMENDMENT-6).
4. The G3 baseline is log(seq_len) + WT log-likelihood. It has no taxon term and uses WTLL, not pseudo-perplexity.
5. The external split is by year (>= 2022), a proxy for the ESM-2 training cutoff.
6. Six assays with empty UniProt annotation (stub entries) were dropped after exclusions.tsv was written (AMENDMENT-5).

## Protocol items NOT delivered (not run)
- ClinVar secondary exposure: not run.
- Mediator check: not run (no code exists).
- AMENDMENT-4 recency sensitivity check (UniProt citations postdating the DMS paper): not run. It could push the external estimate up; this is untested.

## Other caveats
- Scoring ran twice because of an operator error (a re-executed Colab cell restarted it). The run is deterministic; the results come from the complete second run.
- analysis.py was tested only on synthetic random data before any real result existed. One column-name bug was fixed in that test, before outcomes.
