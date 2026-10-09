# DOC-2-073 Protocol (R1 of DOC-2-009) - FROZEN BEFORE ANY OUTCOME

Locked 2026-10-09 IST. No DMS scores, labels or ESM-2 outputs have been read for this project at lock time.
Amendments must be new dated files (AMENDMENT-n.md) made before outcome access and never overwrite this file.

## Question
Is ESM-2 zero-shot variant-effect reliability lower for proteins with thin literature/disease annotation ("understudied") than for heavily studied ones, after controlling for protein length, taxon and protein family?

## Data (open, hash-frozen in DATA_MANIFEST.tsv before labels are read)
- Reliability labels: ProteinGym substitution DMS assays (latest release at acquisition; version + sha256 recorded).
- Annotation depth (exposure): per UniProt entry, log(1 + number of distinct PubMed references) and UniProt annotation score, from one frozen UniProt release (version recorded). ClinVar variant count per gene is a secondary exposure.
- Protein family: InterPro connected components (as in DOC-2-009).
- Model: ESM-2 650M, masked-marginal scoring, fp16 on Colab T4. Proteins >1022 aa use the fixed windowing rule in code, frozen at commit time.

## Outcome per assay
Reliability = Spearman rho between ESM-2 masked-marginal score and DMS fitness score.

## Splits (locked)
- EXTERNAL set: assays whose source publication date is 2022-01-01 or later. Model-card training cutoff to be checked and recorded before labels are read; if ESM-2 training data postdates 2022-01, the cutoff moves later, by amendment, before any outcome. External set is touched once.
- DISCOVERY set: all earlier assays. Cross-validation groups: leave-one-taxon-group-out (human, other eukaryote, prokaryote, virus) and InterPro-component-grouped folds.
- Virus stratum is exploratory only and cannot rescue a failed gate.

## Gates (all must pass for a POSITIVE label)
- G1 (discovery association): partial Spearman between annotation depth and reliability, controlling log length and taxon, >= +0.20, with 95% family-cluster bootstrap CI excluding 0 (10,000 resamples).
- G2 (frozen external): same sign as G1, partial rho >= +0.15, CI excluding 0 - only evaluated if external N >= 40 assays. If N < 40 the label is INCONCLUSIVE-UNDERPOWERED, not positive.
- G3 (predictive utility): adding annotation-depth features to baseline features (log length, taxon, WT pseudo-perplexity) improves leave-taxon-out R2 by >= 0.03 over baseline, and a family-block permutation test gives p < 0.05.
- Direction is pre-declared: more annotation, higher reliability. The opposite sign at any gate is a reported negative, not a new hypothesis.

## Failure policy
- Any failed gate gives a labeled HONEST NEGATIVE with full tables. No gate relaxation, no feature fishing, no dropping assays after seeing results, no swapped datasets.
- Excluded assays are limited to pre-declared reasons (sequence not mappable to UniProt, <50 scored variants), listed before scoring.
- Mediator check (exploratory, labeled as such): whether WT pseudo-perplexity explains any annotation effect.

## Deliverables
Report, per-assay table, code with acquisition commands (no temp paths - fixing the DOC-2-009 reproducibility gap), hashes, and a results README to this repo and the Drive folder.
