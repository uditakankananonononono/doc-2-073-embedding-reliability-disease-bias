# DOC-2-073 results (frozen analysis, AMENDMENT-6; no post-hoc changes)

Label: HONEST NEGATIVE. No gate passed. No claim of an annotation-bias effect is made.

- Scoring: ESM-2 650M wt-marginal, fp32, Colab CPU, 204 assays, files per_assay.csv (rho per assay, WT log-likelihood, length).
- G1 (discovery, year <= 2021, n=98 assays after n_scored >= 50 filter): partial rho = -0.068, 95% protein-cluster bootstrap CI [-0.251, 0.140]. Fail.
- G2 (external, year >= 2022, n=106): partial rho = +0.101, CI [-0.137, 0.318]. Fail (also opposite sign to G1; CIs span 0).
- G3 (leave-one-taxon-out on Human/Eukaryote/Prokaryote; baseline log(seq_len)+WTLL): baseline R2 -0.210, full R2 -0.252, delta -0.042, permutation p = 0.544. Fail. Virus fold (exploratory only): base 0.280, full 0.148.
- Caveats: UniProt_ID is the cluster unit (InterPro clustering not acquired, AMENDMENT-6). Scoring was run twice because of an operator error (second run overwrote the first; deterministic, results used are from the complete second run). analysis.py was tested only on synthetic random data before any real result existed.
- Independent review of these results is pending before any claim.
