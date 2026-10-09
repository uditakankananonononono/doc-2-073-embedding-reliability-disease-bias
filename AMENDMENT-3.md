# AMENDMENT-3 (2026-10-09 IST) - pre-outcome clarifications requested by the independent gate

1. WTLL is the mean over ALL residues of the full target sequence (the ProteinGym target_seq, or the mapped full UniProt sequence if used), computed once per protein and fixed per protein, never per assay and never over only the scored positions.
2. Multi-mutant variants: excluded from scoring. Only single amino-acid substitutions are scored; no additive extension is used. The "fewer than 50 scored variants" exclusion counts single substitutions only, so an assay with fewer than 50 single substitutions is excluded, and the exclusion is recorded in the pre-declared exclusion list before scoring.
3. Runtime: the RSS measurement in DATA_MANIFEST.tsv is taken in the same Colab CPU runtime session that performs the scoring (runtime type, total RAM and peak RSS recorded). If the scoring session is a different runtime, the measurement is repeated in it before labels are read and the model rule is re-applied.
Everything else in PROTOCOL.md, AMENDMENT-1 and AMENDMENT-2 stands.
