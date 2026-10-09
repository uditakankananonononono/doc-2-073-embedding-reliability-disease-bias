# AMENDMENT-4 (2026-10-09 IST) - pre-outcome data-definition notes, written before any label is read

1. External split: ProteinGym's reference file gives publication `year` only. "Source publication date >= 2022-01-01" is implemented as year >= 2022. Discovery = year <= 2021.
2. "Sequence not mappable to UniProt" is implemented as: the ProteinGym UniProt_ID entry is not found in the frozen UniProt query. No sequence-identity test is used for exclusion.
3. Exclusions are decided from reference-file metadata only (single-substitution count < 50; entry not found) and listed in exclusions.tsv before any DMS score file is opened.
4. Exposures: n_pubmed = distinct PubMed IDs cited in the UniProt entry (lit_pubmed_id); annotation_score = UniProt annotation score. ClinVar counts are NOT acquired; no gate depends on them.
5. Known limitation, stated now: UniProt citations are counted at the frozen release date, which can include papers published after a DMS assay (including the DMS paper itself). This can inflate n_pubmed for recent DMS targets. It is reported as a limitation and tested only as an exploratory sensitivity check; it does not change any gate.
6. acquire_data.py reads only metadata and never opens a DMS score file. Hashes go to DATA_HASHES.tsv.
